import os, re
from ij import IJ, WindowManager, Prefs
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import Roi, Toolbar, WaitForUserDialog
import os
from ij.io import DirectoryChooser, OpenDialog
from ij.process import ImageProcessor
from ij.plugin import ChannelSplitter

# Bioformats
from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions

# get input path for merge file
opener = DirectoryChooser("Select the input folder")
input_folder = opener.getDirectory()
print input_folder


#Define results output folder
folder = input_folder+"timepoints/"
if not os.path.exists(folder):
	os.mkdir(folder)

# collect list of files to be processed
#file_list = [filename for filename in os.listdir(input_folder) if ".tif" in filename]
# Recursive option
file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(input_folder) for f in fn if '.tif' in f]

for file_path in file_list:

	# read in and display ImagePlus(es) with arguments
	options = ImporterOptions()
	options.setId(file_path)
	options.setSplitTimepoints(True)

	imps = BF.openImagePlus(options)
	for imp in imps:
	    imp.show()
	    name = imp.getTitle()
	    details = name.split("_")[0]
	    print file_path
	    pos = details.split('Pos')[1]
	    time = name.split("=")[1]
	    image_type = re.split(r'.lif |- |_|\\', file_path)[-3]
	    mutant = re.split(r'.lif |- |_|\\', file_path)[-6]
	    cotrans = re.split(r'.lif |- |_|\\', file_path)[-5]
	    new_name = mutant+"_"+cotrans+"_pos_"+pos+"_t_"+time+"_"+image_type
	    print new_name
	    IJ.saveAs(imp, "Tiff", os.path.join(folder, new_name))
	    imp.close()

	print "Done!"
