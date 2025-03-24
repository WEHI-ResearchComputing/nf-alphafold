""" Proetein Fold Utilities class """

import base64
import logging
import glob
import json
import math
import os
import pickle
import re
import numpy as np
import matplotlib.pyplot as plt # type: ignore
from typing import Tuple

MODEL_PRESETS = {
    'monomer': (
        'model_1',
        'model_2',
        'model_3',
        'model_4',
        'model_5',
    ),
    'monomer_ptm': (
        'model_1_ptm',
        'model_2_ptm',
        'model_3_ptm',
        'model_4_ptm',
        'model_5_ptm',
    ),
    'multimer': (
        'model_1_multimer_v3',
        'model_2_multimer_v3',
        'model_3_multimer_v3',
        'model_4_multimer_v3',
        'model_5_multimer_v3',
    ),
}
MODEL_PRESETS['monomer_casp14'] = MODEL_PRESETS['monomer']

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class ProteinFold:
    """ ProteinFold Class defined to generate plots, html reports and rank results"""
    def __init__(self, sample_name:str, fold_type:str, structures:str,
                output_dir:str,
                model_preset:str,
                num_predictions_per_model:int):
        self.sample_name = sample_name
        self.fold_type = fold_type
        self.structures=structures
        self.output_dir=output_dir
        self.num_predictions_per_model=num_predictions_per_model
        self.model_preset=model_preset
        self.model_names_extra = [f'{m}_pred_{i}' for m in MODEL_PRESETS[self.model_preset] for i in range(self.num_predictions_per_model)]
        self.plots_dir = os.path.join(self.output_dir, "plots")

    def generate_html_from_template(self, html_template:str, mean_plddt:float):
        """Function generating the final html report from the template"""
        logging.info("Generating html report...")

        with open(html_template, "r",encoding="utf-8") as htmlfile:
            template=htmlfile.read()
            template = template.replace("*sample_name*", self.sample_name )
            template   = template.replace("*prog_name*", self.fold_type)
            args_pdb_array_js = ",\n".join([f'"{os.path.basename(model)}"' for model in self.structures])
            template = re.sub(
                r'const MODELS = \[.*?\];',  # Match the existing MODELS array in HTML template
                f'const MODELS = [\n  {args_pdb_array_js}\n];',  # Replace with the new array
                template,
                flags=re.DOTALL,
            )

            averages_plddt = 'const LDDT_AVERAGES = [{}];'.format(",\n".join(map(str,mean_plddt)))
            template = template.replace("const LDDT_AVERAGES = [];", averages_plddt)

            for i, structure in enumerate(self.structures):
                with open(structure, "r", encoding="utf-8") as struct_file:
                    template = template.replace(f"ranked_{i}.pdb*", struct_file.read().replace("\n", "\\n"))
            self._embed_images_in_template(template)
            logging.info("HTML report generated successfully.")


    def _embed_images_in_template(self, template:str):
        """Embed images into the HTML template"""
        image_files = [
            ("seq_coverage.png", f"{self.output_dir}/plots/{self.sample_name}_coverage_LDDT.png"),
            ("pae.png", f"{self.output_dir}/plots/{self.sample_name}_Pae.png"),
            ("plddt.png", f"{self.output_dir}/plots/{self.sample_name}_PLDDT.png")
        ]

        for placeholder, image_path in image_files:
            with open(image_path, "rb") as in_file:
                encoded_image = base64.b64encode(in_file.read()).decode('utf-8')
                template = template.replace(placeholder, f"data:image/png;base64,{encoded_image}")

        with open(f"{self.output_dir}/{self.sample_name}_report.html", "w", encoding="utf-8") as out_file:
            out_file.write(template)
        logging.info("Images embedded successfully.")

    def calculate_plddts(self) -> dict:
        """Calculate pLDDT scores for each model"""
        plddts = {}
        for model_name in self.model_names_extra:
            for filenm in glob.glob(os.path.join(self.output_dir, f'result_{model_name}.pkl')):
                logging.info("Found results for %s", model_name)
                with open(filenm, 'rb') as result_file:
                    result = pickle.load(result_file)
                    plddts[model_name] = np.mean(result['plddt'])

        return plddts

    def rank_models(self, plddts: dict) -> list:
        """Rank models based on pLDDT scores"""
        ranked_order = []
        for idx, (model_name, _) in enumerate(sorted(plddts.items(), key=lambda x: x[1], reverse=True)):
            ranked_order.append(model_name)
            ranked_output_path = os.path.join(self.output_dir, f'ranked_{idx}.pdb')
            with open(ranked_output_path, 'w', encoding='utf-8') as f:
                relaxed_pdb_path = os.path.join(self.output_dir, f'relaxed_{model_name}.pdb')
                unrelaxed_pdb_path = os.path.join(self.output_dir, f'unrelaxed_{model_name}.pdb')
                if os.path.exists(relaxed_pdb_path):
                    with open(relaxed_pdb_path, 'r', encoding='utf-8') as relaxed_pdbs:
                        f.write(relaxed_pdbs.read())
                        logging.info("For %s - relaxed model saved", model_name)
                elif os.path.exists(unrelaxed_pdb_path):
                    with open(os.path.join(self.output_dir,f'unrelaxed_{model_name}.pdb'),'r', encoding='utf-8') as unrelaxed_pdbs:
                        f.write(unrelaxed_pdbs.read())
                        logging.info("For %s - unrelaxed saved", model_name)
        ranking_output_path = os.path.join(self.output_dir, 'ranking_debug.json')
        with open(ranking_output_path, 'a', encoding='utf-8') as f:
            f.write(json.dumps({'plddts': plddts, 'order': ranked_order}, indent=4))
        
        return ranked_order

    def plot_results(self, ranked_order: list, plddts: dict):
        """Plot pLDDT and PAE results"""
        os.makedirs(self.plots_dir, exist_ok=True)
        feature_dict = pickle.load(open(os.path.join(self.output_dir, "features.pkl"),'rb'))
        self._plot_seq_coverage(feature_dict)
        self._plot_plddt(ranked_order, plddts)
        self._plot_pae(ranked_order, plddts)

    def _plot_seq_coverage(self, feature_dict: dict):
        msa = feature_dict['msa']
        seqid = (np.array(msa[0] == msa).mean(-1))
        seqid_sort = seqid.argsort()
        non_gaps = (msa != 21).astype(float)
        non_gaps[non_gaps == 0] = np.nan
        final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
        plt.figure(figsize=(14, 4), dpi=100)
        plt.imshow(final,
                interpolation='nearest', aspect='auto',
                cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
        plt.plot((msa != 21).sum(0), color='black')
        plt.xlim(-0.5, msa.shape[1] - 0.5)
        plt.ylim(-0.5, msa.shape[0] - 0.5)
        plt.colorbar(label="Sequence identity to query", )
        plt.xlabel("Positions")
        plt.ylabel("Sequences")
        plt.savefig(f"{self.output_dir}/plots/{self.sample_name}_coverage_LDDT.png")

    def _plot_plddt(self, ranked_order: list, plddts: dict):
        ranking_dict={'plddts': plddts, 'order': ranked_order}
        model_dicts={}
        for model_name in self.model_names_extra:
            for filenm in glob.glob(os.path.join(self.output_dir,f'result_{model_name}.pkl')):
                result=pickle.load(open(filenm, 'rb'))
                model_dicts[model_name] = result

        s = 0
        plt.figure(figsize=(14, 4), dpi=100)
        for model_name, value in model_dicts.items():
            plt.plot(value["plddt"], 
                    label=f"{model_name}, plddts: {round(list(ranking_dict['plddts'].values())[s], 6)}")
            s += 1
        #plt.legend()
        plt.legend(loc='lower left')
        plt.ylim(0, 100)
        plt.ylabel("Predicted LDDT")
        plt.xlabel("Positions")
        plt.savefig(f"{self.output_dir}/plots/{self.sample_name}_PLDDT.png")

    def _plot_pae(self, ranked_order: list, plddts: dict):
        ranking_dict={'plddts': plddts, 'order': ranked_order}
        model_dicts={}
        for model_name in self.model_names_extra:
            for filenm in glob.glob(os.path.join(self.output_dir,f'result_{model_name}.pkl')):
                result=pickle.load(open(filenm, 'rb'))
                model_dicts[model_name] = result
        plt.figure(figsize=(14, 4), dpi=100)
        if "predicted_aligned_error" in model_dicts[self.model_names_extra[0]]:
            for n, (model_name, value) in enumerate(model_dicts.items()):
                plt.figure(figsize=[8 * 2, 6],dpi=100)
                plt.subplot(1, 2, 1)
                plt.title(f'Predicted LDDT, {model_name}')
                plt.plot(value["plddt"], 
                label=f"{model_name}, plddts: {round(list(ranking_dict['plddts'].values())[n], 6)}")
                plt.subplot(1, 2, 2)
                plt.title(f'Predicted Aligned Error, {model_name}')
                plt.imshow(value["predicted_aligned_error"], label=model_name, vmin=0., vmax=value["max_predicted_aligned_error"], cmap='Greens_r')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.savefig(f"{self.output_dir}/plots/{self.sample_name}_Pae.png")
        else:
            logging.info("No predicted_aligned_error found. Make sure you choose monomer_ptm when running AlphaFold prediction.")