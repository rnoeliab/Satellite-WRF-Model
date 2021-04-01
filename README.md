# Satellite-WRF-Model
This project is create to work with satellite and model data. The programming language is python and the operating system is linux.

## MODIS
### How to get the MODIS data?

* To start using MODIS data, you need to go to the [MODIS data](https://ladsweb.modaps.eosdis.nasa.gov/search/) website. But first you need to register. So, let's go to the following link: [REGISTER](https://urs.earthdata.nasa.gov/users/new?client_id=A6th7HB-3EBoO7iOCiCLlA&redirect_uri=https%3A%2F%2Fladsweb.modaps.eosdis.nasa.gov%2Fcallback&response_type=code&state=L3NlYXJjaC8%3D%0A).
* When you are already registered , you will find five steps to follow: 
```
  1. PRODUCTS --> you need to select one sensor (MODIS) and satellite (Terra or Aqua), considering the collection 6.1. Then the place of study (in this case, "Atmosphere") and the variable (Aerosol). With this information, you will see two options (MOD04_3K\MYD04_3K  and MOD04_L2\MYD04_L2). The difference is because "3K" is to 3km spatial resolution and "L2" is to 10km spatial resolution. Of these two options you can select one or two or more if you have the option. 
  2. TIME --> you need to select the period of study and click in "add Date".
  3. LOCATION --> In this part, select area of interest by "Tiles". For me this option is easier. 
  4. FILES --> Here, you will see a list of all the images that contain the product that was chosen, within the period and the study area. So, you can download one by one (clicking in the last columna) or select various (clicking on the name of the images)or all images.
  If you are registered and you did "login" you will go to the last step, but if you did not "login" you will only have to confirm your username and password. 
  5. REVIEW & ORDER --> Finally, you must confirm the order and click on "Submit Order".
```
* The order of the images takes mostly minutes to arrive in the registered e-mail. 
* When you get the email "LAADS Web Order Notification", there you will find the order number: Your Order ID is: 501481680
* When you have access to the order file, the following step is download the MODIS images. For this, I am going to help you do that by leaving you the following steps:
```
1. Read the ["download_MODIS.sh"] script, in this script, find the path where the file file_order.txt could be found, if it doesn't exist create it.
         Ex: /home/username/file_order.txt
In the file file_order.txt, write "/archive/orders/501481680/" (without the quotes), where "501481680" is the order number that arrived in your email.
2. On the MODIS data website, where you performed the five steps to download the MODIS images, in the upper right, in "Profile", click on "App Keys" (assuming you are registered and did login), on "App Keys" you will find your key for MODIS. 
        Ex: "10E74C60-F16A-11E8-B781-B0B1502E57AA"
copy your key and paste it on line 7 of the "download_MODIS.sh" script, after the "Bearer" command.
3. Finally, put the path where you want to download the images. 
```
The download_MODIS script is located here:[download_MODIS.sh](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/download_MODIS.sh)
```shell
#! /bin/bash
# read the order files 
filename='/home/username/file_order.txt'

while read line; do 
echo $line
wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov$line" --header "Authorization: Bearer 10E74C60-F16A-11E8-B781-B0B1502E57AA" -P /home/username/output/
done< $filename 

echo "Finished the download"
``` 
* To run this script, do the following: In the terminal, write
```
chmod +x download_MODIS.sh
./download_MODIS.sh
```


## OMI
### How to get the MODIS data?
In this part, we can start downloading the OMI data using the following link: [GES DISC](https://disc.gsfc.nasa.gov/datasets/).
