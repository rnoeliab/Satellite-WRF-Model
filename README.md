# Satellite-WRF-Model
This project is create to work with satellite and model data. The programming language is python.

## MODIS
### How to get the MODIS data?
* To start using MODIS data, you need to go to the [MODIS data](https://ladsweb.modaps.eosdis.nasa.gov/search/) website. But first you need to register. So, let's go to the following link: [REGISTER](https://urs.earthdata.nasa.gov/users/new?client_id=A6th7HB-3EBoO7iOCiCLlA&redirect_uri=https%3A%2F%2Fladsweb.modaps.eosdis.nasa.gov%2Fcallback&response_type=code&state=L3NlYXJjaC8%3D%0A).
* When you are already registered , you will find five steps to follow: 
  PRODUCTS --> you need to select one sensor (MODIS) and satellite (Terra or Aqua), considering the collection 6.1. Then the place of study (in this case, "Atmosphere") and the variable (Aerosol). With this information, you will see two options (MOD04_3K\MYD04_3K  and MOD04_L2\MYD04_L2). The difference is because "3K" is to 3km spatial resolution and "L2" is to 10km spatial resolution. Of these two options you can select one or two or more if you have the option. 
  TIME --> you need to select the period of study and click in "add Date".
  LOCATION --> In this part, select area of interest by "Tiles". For me this option is easier. 
  FILES --> Here, you will see a list of all the images that contain the product that was chosen, within the period and the study area. So, you can download one by one (clicking in the last columna) or select various (clicking on the name of the images)or all images.
  If you are registered and you did "login" you will go to the last step, but if you did not "login" you will only have to confirm your username and password. 
  REVIEW & ORDER --> Finally, you must confirm the order and click on "Submit Order".
* The order of the images takes mostly minutes to arrive in the e-mail. 


## OMI
