import zipfile
import os
import shutil

base_dir = r"d:\Cache\Yandex.Disk\РАЗРАБОТКА\code\imerg2meteo_calib\data"
meteo_zip = os.path.join(base_dir, "срочные данные_осадки.zip")
imerg_zip = os.path.join(base_dir, "IMERG_STATIONS_3H.zip")

meteo_out = os.path.join(base_dir, "meteo")
imerg_out = os.path.join(base_dir, "imerg")

if os.path.exists(meteo_out):
    shutil.rmtree(meteo_out)
if os.path.exists(imerg_out):
    shutil.rmtree(imerg_out)

os.makedirs(meteo_out, exist_ok=True)
os.makedirs(imerg_out, exist_ok=True)

if os.path.exists(meteo_zip):
    print("Extracting meteo data with cp866 encoding...")
    with zipfile.ZipFile(meteo_zip, 'r', metadata_encoding='cp866') as z1:
        z1.extractall(meteo_out)

if os.path.exists(imerg_zip):
    print("Extracting IMERG data with cp866 encoding...")
    with zipfile.ZipFile(imerg_zip, 'r', metadata_encoding='cp866') as z2:
        z2.extractall(imerg_out)

print("Extraction complete.")
