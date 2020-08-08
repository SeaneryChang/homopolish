import argparse
from version import __version__
from modules.polish_interface import polish_genome
from modules.train_interface import train_model

def add_polish_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-a",
        "--assembly",
        type=str,
        required=True,
        help="[REQUIRED] Path to a assembly genome."
    )
    parser.add_argument(
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="[REQUIRED] Path to a trained model (pkl file). Please see our github page to see options."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-s",
        "--sketch_path",
        type=str,
        required=False,
        help="Path to a mash sketch file."
    )
    group.add_argument(
        "-g",
        "--genus",
        type=str,
        required=False,
        help="Genus name"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use. [1]"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=False,
        default='./output/',
        help="Path to the output directory. [output]"
    )
    parser.add_argument(
        "--minimap_args",
        type=str,
        required=False,
        default='asm5',
        help="Minimap2 -x argument. [asm5]"
    )
    parser.add_argument(
        "--mash_threshold",
        type=str,
        required=False,
        default='0.95',
        help="Mash output threshold. [0.95]"
    )
    parser.add_argument(
        "--download_contig_nums",
        type=str,
        required=False,
        default='20',
        help="How much contig to download from NCBI. [20]"
    )
    parser.add_argument(
        "-d",
        "--debug",
        required=False,
        action = "store_false",
        help="Keep the information of every contig after mash, such as homologous sequences and its identity infomation. [no]"
    )
    parser.add_argument(
        "--mash_screen",
        required=False,
        action='store_true',
        default=False,
        help="Use mash screen. [mash dist]"
    )

    return parser
def add_train_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-d",
        "--dataframe_dir",
        type=str,
        required=True,
        help="[REQUIRED] Path to a directory for alignment dataframe."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=False,
        default='./output/',
        help="Path to the output directory. [output]"
    )
    parser.add_argument(
        "-p",
        "--output_prefix",
        type=str,
        required=False,
        default='train',
        help="Prefix for the train model. [train]"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use. [1]"
    )
    return parser

def main():
    parser = argparse.ArgumentParser(description="Homopolish is a SVM based polisher for polishing ONT-based assemblies. \n"
                                                 "1) polish: Run the polishing pipeline.\n"
                                                 "2) train: Train your own SVM model.",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    subparsers = parser.add_subparsers(dest='sub_command')
    parser_polish = subparsers.add_parser('polish', help="Run the polishing pipeline.")
    add_polish_arguments(parser_polish)

    parser_train = subparsers.add_parser('train', help="Train a model.")
    add_train_arguments(parser_train)

    FLAGS, unparsed = parser.parse_known_args()
    if FLAGS.sub_command == 'polish':        
        polish_genome(FLAGS.mash_screen, FLAGS.assembly, FLAGS.model_path, FLAGS.sketch_path, FLAGS.genus, FLAGS.threads, \
                FLAGS.output_dir, FLAGS.minimap_args, FLAGS.mash_threshold, FLAGS.download_contig_nums, FLAGS.debug)

    elif FLAGS.sub_command == 'train':
        train_model(FLAGS.dataframe_dir, FLAGS.output_dir, FLAGS.output_prefix, FLAGS.threads)

    elif FLAGS.version is True:
        print("Homopolish VERSION: ", __version__)

    
if __name__ == "__main__":
    main()