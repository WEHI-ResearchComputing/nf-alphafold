#!/stornext/System/data/apps/rc-tools/rc-tools-1.0/bin/tools/envs/biopython/bin/python
"""
 Script to generate plots and interactive visualisaion 
 of AlphaFold results.
"""
import argparse
import os
from collections import OrderedDict
import logging
import matplotlib.pyplot as plt
import numpy as np

from proteinfold import ProteinFold

##import plotly.graph_objects as go


logging.disable(logging.DEBUG)

def main():
    """Main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb',   dest='pdb',required=True, nargs="+")
    parser.add_argument('--name',  dest='name')
    parser.add_argument('--output_dir',dest='output_dir')
    parser.add_argument('--html_template',dest='html_template')
    parser.add_argument('--model_preset',dest='model_preset')
    parser.set_defaults(output_dir='./')
    parser.set_defaults(html_template='../assets/proteinfold_template.html')
    parser.set_defaults(name='sample')
    args = parser.parse_args()

    structures = args.pdb
    structures.sort()
    fold=ProteinFold(sample_name=args.name,
                    fold_type="AlphaFold",
                    structures=structures,
                    output_dir=args.output_dir,
                    num_predictions_per_model=1,
                    model_preset=args.model_preset)
    
    plddts = fold.calculate_plddts()
    ranked_order = fold.rank_models(plddts)
    logging.info("Ranked models: %s", ranked_order)
    #_, max_plddt= (max(sorted(plddts.items(), key=lambda x: x[1], reverse=True)))
    mean_plddt = list(plddts.values())
    #max_plddt = max(plddts.values())
    logging.info("Mean pLDDT: %s", mean_plddt)
    fold.plot_results(ranked_order, plddts)
    fold.generate_html_from_template(html_template=args.html_template,mean_plddt=mean_plddt)
    

if __name__ == "__main__":
    main()