README for Term_project
Folder structure:
Code
mtbs_processing
mtbs_combined
mtbs_unzipped
temp
mtbs_utm
mtbs_ufi_combined
raw_data
burn_severity_data
mtbs_data
ufi_data
response_variable_data
uwrap_data
wrc_data

Project objective:
Burn severity is commonly assessed in terms of the magnitude of ecological change post-fire. Burn severity is an important indicator of landscape response to wildfire and there exist programs and agencies which aim to quantify and catalogue post-fire burn severity. Post-fire burn severity data is commonly available in the form of raster layers containing data on the Relative difference Normalized Burn Ratio (RdNBR), or the difference between pre- and post-fire NBR. These data are available from public repositories, such as the Monitoring Trends in Burn Severity (MTBS) databases, or private researches, such as the Utah Fire Atlas (UFI), assembled by Megan Nasto and Jim Lutz out of Utah State University. 
While these data are publicly available for large wildfires (>500 acres) from MTBS, the file structure limits the accessibility of assembling large MTBS datasets, as data for each fire is contained in a unique zipped folder. The goal of this project is to create a file directory to organize intermediate processing files and final products for extracting MTBS burn severity data for wildfires in the state of Utah between 2012-2018. The specific objective of this project is to produce a raster layer of RdNBR burn severity data for wildfires in the state of Utah from 2012-2018. This data will be used to assess the efficacy of decision support tool (DST) data which offer spatial information on fire risk or threat such as the Utah Wildfire Risk Assessment Portal (UWRAP) or the Wildfire Risk to Communities (WRC) website. 
While MTBS data require intermediate processing (unzipping, extracting, projecting, etc.), UFI data, however, came with a combined layer of smaller (500>100 acres) wildfires in the state of Utah. This will require less intermediate work before being combined with MTBS data layers. I have chosen a project-based directory structure for this project as it is self-contained. The final products of this processing will be stored in the folder ‘mtbs_ufi_combined’ where the clipped, projected and combined raster layers from mtbs and ufi will be combined into a discrete layer of RdNBR burn severity data for wildfires >100 acres for the state of Utah from 2012-2018. 
