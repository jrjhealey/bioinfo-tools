# Extract beginning and end co-ordinates from the output
# of MultiGeneBlast. This script will produce a table of
# starts and ends, that can be used with the genbank slicer
# script to subset operons for each resulting hit.


import traceback, warnings

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "MGBparser"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"

def get_args():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Extract hit beginnings and ends from the output of MultiGeneBlast.')

        parser.add_argument('clusterfile',
                            action='store',
                            help='The text file of hits output by MGB. By default this is called \'clusterblast_output.txt\'.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()



def main():
	args = get_args()






if __name__ == "__main__":
	main()
