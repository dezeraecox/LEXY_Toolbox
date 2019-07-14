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

def close_images():
	open_images = WindowManager.getImageTitles()
	for imagename in open_images:
		print imagename
		IJ.selectWindow(imagename)
		IJ.run("Close")
	print "Windows closed"

def nuclei_processor(imp_particles, thresh_type, folder, impname, channel_name):
	rm = RoiManager.getInstance()
	if not rm:
  		rm = RoiManager()
	rm.reset()

	# define a results table
	rt = ResultsTable.getResultsTable()

	imp_particles.show()
	# generate thresholded image for ROI
	nuclei = imp_particles.duplicate()
	nuclei.show()
	IJ.run("Gaussian Blur...", "sigma=3")
	IJ.setAutoThreshold(nuclei, thresh_type)
	IJ.run("Convert to Mask")
	IJ.run("Fill Holes")
	IJ.run("Watershed")

	# select thresholded image (to set ROIs)
	IJ.run("Set Measurements...", "area mean standard min area_fraction limit display add redirect=["+imp_particles.title+"] decimal=3")
	IJ.run("Analyze Particles...", "size=30-Infinity show=Outlines display clear add")

	# get the ROI manager and save
	rm.runCommand("save selected", os.path.join(folder, impname+'_'+channel_name+"_ROIs.zip"))

#	pixel_collector(rm, imp_measure, channel_name, impname, folder)
#	return rm

def cell_processor(imp_particles, imp_measure, folder, impname, channel_name):

	# Reset ROI manager
	rm = RoiManager.getInstance()
	if not rm:
			rm = RoiManager()
	rm.reset()

	# select image to measure, duplicate and show
	imp_measure.show()
#	IJ.run("Duplicate...", "title=Overlay_measure")
#	IJ.run("Red")

	imp_particles.show()
	IJ.run("Duplicate...", "title=DIC")
#	IJ.run("Add Image...", "image=[Overlay_measure] x=0 y=0 opacity=20")

	pause = WaitForUserDialog("Select centre of all displayed cells \n \nPress OK to continue")
	pause.show()
	IJ.run("Cell Outliner", "cell_radius=28 tolerance=1.0 kernel_width=13 dark_edge kernel_smoothing=1 polygon_smoothing=1 weighting_gamma=3 iterations=3 dilate=0")

	dic_outline = WindowManager.getImage("DIC Cell Outline")
	#WindowManager.setWindow(dic_outline)
	dic_outline.getProcessor().setAutoThreshold("Default dark")
	IJ.setThreshold(1, 255, "Black & White")
	IJ.run("Convert to Mask")
	IJ.run("Watershed")

	print "Measuring cells in "+imp_measure.title
	IJ.run("Set Measurements...", "area mean standard min area_fraction limit display add redirect=["+imp_measure.title+"] decimal=3")
	IJ.run("Analyze Particles...", "size=30-Infinity show=Outlines display clear add")

	# get the results table and save
	rt = ResultsTable.getResultsTable()
	rt.saveAs(os.path.join(folder, impname+"_cells_Results.csv"))
	print "Analyse particles completed for cells"

	# get the ROI manager and save
	rm = RoiManager.getInstance()
	rm.runCommand("save selected", os.path.join(folder, impname+"_cells_ROIs.zip"))

	pixel_collector(rm, imp_measure, channel_name, impname, folder)

	return rm


def pixel_collector(rm, channel_imp, channel_name, impname, folder):

		# define new Results table
		rt = ResultsTable()

		IndRois = rm.getIndexes()
		for index in IndRois:
			ROI = rm.getRoi(index)
			ROI_name = ROI.getName()
			coords = ROI.getContainedPoints()

			row = 0
			for pixel in coords:
				x_coord = pixel.getX()
				y_coord = pixel.getY()

				rt.setValue(ROI_name+"_X_pos", row, int(x_coord))
				rt.setValue(ROI_name+"_Y_pos", row, int(y_coord))

				pixel_2 = channel_imp.getProcessor().getPixel(int(x_coord), int(y_coord))
				rt.setValue(ROI_name+"_"+channel_name, row, pixel_2)

		  		row = row + 1
		rt.show("Results")

		rt.save(os.path.join(folder, impname+'_'+channel_name+"_pixels.csv"))
		print "Pixel collection done!"

## ------------- processing happens here -----------------
# get input path for merge file
opener = DirectoryChooser("Select the input folder")
input_folder = opener.getDirectory()

#Define results output folder
folder = input_folder+"Results\\"
if not os.path.exists(folder):
	os.mkdir(folder)
# collect list of files to be processed
file_list = [filename for filename in os.listdir(input_folder) if ".tif" in filename]

print file_list

before_list = [file_name for file_name in file_list if 'Before' in file_name]
after_list = [file_name for file_name in file_list if 'Activation' in file_name]
t0_list = [file_name for file_name in after_list if 't_0' in file_name]

# Define ROI manager, clear if already exists
rm = RoiManager.getInstance()
if not rm:
		rm = RoiManager()
rm.reset()


# Process cells in t0_activation images
for impfile in t0_list:
	# open image
	imp = IJ.openImage(os.path.join(input_folder, impfile))
	imp.show()

	impname=impfile.strip('.tif')
	imp_details = impname.split('_')
	mutant, cotrans, _, pos_number, _, time_point, image_type = imp_details

	print mutant, pos_number, time_point


	imps = ChannelSplitter.split(imp)
	cell_processor(imp_particles=imps[2], imp_measure=imps[1], folder=folder, impname=impname, channel_name='mCherry_cells')

	# Process cell ROIs
	cell_ROIs = impname + '_cells_ROIs.zip'
	print cell_ROIs
	rm.runCommand("Open", os.path.join(input_folder, 'Results', cell_ROIs))
	pixel_collector(rm, imps[4], 'GFP_cells', impname, folder)

	close_images()



# Process after files, collecting cell and nuclei info
for impfile in after_list:
	# open image
	imp = IJ.openImage(os.path.join(input_folder, impfile))
	imp.show()

	impname=impfile.strip('.tif')
	imp_details = impname.split('_')
	mutant, cotrans, _, pos_number, _, time_point, image_type = imp_details

	print mutant, pos_number, time_point

	## Process nuclei image
	imps_thresh = ChannelSplitter.split(imp)
	nuclei_processor(imp_particles=imps_thresh[3], thresh_type='Li dark', folder=folder, impname=impname, channel_name='nuclei')

	imps = ChannelSplitter.split(imp)
	# Process nuc ROIs in mCherry
	rm.reset()
	nuc_ROIs = impname + '_nuclei_ROIs.zip'
	print nuc_ROIs
	rm.runCommand("Open", os.path.join(input_folder, 'Results', nuc_ROIs))
	pixel_collector(rm, imps[1], 'mCherry_nuclei', impname, folder)
#	close_images()
	
	# Process nuc ROIs in mCherry
	rm.reset()
	nuc_ROIs = impname + '_nuclei_ROIs.zip'
	print nuc_ROIs
	rm.runCommand("Open", os.path.join(input_folder, 'Results', nuc_ROIs))
	pixel_collector(rm, imps[4], 'GFP_nuclei', impname, folder)
	close_images()



# Process before files using the predefined ROI managers to gather cell and nuclei info
for impfile in before_list:
	# open image
	imp = IJ.openImage(os.path.join(input_folder, impfile))
	imp.show()

	impname=impfile.strip('.tif')
	imp_details = impname.split('_')
	mutant, cotrans, _, pos_number, _, time_point, image_type = imp_details

	print mutant, pos_number, time_point

	imps = ChannelSplitter.split(imp)

	# Process cell ROIs
	rm.reset()
	cell_ROIs = impname.strip('_Before') + '_Activation_cells_ROIs.zip'
	print cell_ROIs
	rm.runCommand("Open", os.path.join(input_folder, 'Results', cell_ROIs))
	pixel_collector(rm, imps[0], 'mCherry_cells', impname, folder)

	# Process nuclei ROIs
	rm.reset()
	nuc_ROIs = impname.strip('_Before') + '_Activation_nuclei_ROIs.zip'
	print nuc_ROIs
	rm.runCommand("Open", os.path.join(input_folder, 'Results', nuc_ROIs))
	pixel_collector(rm, imps[0], 'mCherry_nuclei', impname, folder)

	close_images()


print "Process Complete"
print "Results for "+mutant+'_'+pos_number+". Saved in directory: "+folder
