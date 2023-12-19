# image-pairs-stereo-capacity
Estimate the stereoscopic capacity (B/H ratio) of pairs of images from the Pléiades, Pléiades Neo or SPOT6-7 satellites

Very high resolution satellites, such as Pléiades and SPOT6-7, can acquire images in stereo or tri-stereo modes. Still, two monoscopic acquisitions can be combined to produce a digital surface model (DSM). howstereo.py can help you to choose which pair of images offers the best configuration in terms of stereoscopic capacity.

For those agile satellites, the incidence angle, which is used to calculate the B/H ratio, is decomposed in the scan and ortho-scan directions (see the Pléiades and SPOT6-7 user guides for more details). For each acquisition, the scan, ortho-scan and scan azimuth angles are provided in the Airbus archive catalog (http://www.intelligence-airbusds.com/geostore/).

Of course, you should also look at other parameters for choosing the images, such as their common footprint, cloud and snow cover, the time difference between the two acquisitions (which induces changes between the two images and therefore potential correlation issues).

*Thanks to Julien & Charlotte for the precious help!*

### Usage:

```
howstereo.py [-h] [--incid1 ALONG ACROSS] [--az1 AZIMUTH] [--incid2 ALONG ACROSS]
             [--az2 AZIMUTH] [--az_mode {scan,target}] [--input_file FILE] [--show_plot]
```
Optional arguments:
```
-h, --help           show this help message and exit
--inc1 ALONG ACROSS  incidence angles for image 1 (in degrees)
--az1 AZIMUTH        azimuth angle for image 1 (in degrees)
--inc2 ALONG ACROSS  incidence angles for image 2 (in degrees)
--az2 AZIMUTH        azimuth angle for image 2 (in degrees)
--az_mode {scan, target}
                     type of azimuth angle (see the note below; default: scan)
--input_file FILE    input from a file instead (in csv format: inc_scan,inc_ortho,az)
--show_plot          show a 3D, interactive plot
```
Note on azimuth angle:

> TL;DR: if the angle is from a SPOT6|7 DIMAP file, select "target" for --az_mode. Otherwise, select "scan" (default).
> 
> The azimuth of the scan axis (i.e. the angle between geographic north and the image line axis on the ground) is used in the B/H calculation. In the Geostore, this angle is called the Orientation angle. In the DIMAP file, for Pléiades, it corresponds to the AZIMUTH_ANGLE, but for SPOT6|7, the AZIMUTH_ANGLE refers to the target azimuth. For Pléiades Neo, it corresponds to the IMAGE_ORIENTATION. If you provide the target azimuth as input, you must select "target" for --az_mode, and the program will perform the conversion. Otherwise, select "scan".

### Examples:

#### Example 1:

- Image 1 has an incidence of -5.44° in the scan direction, 7.07° in the ortho-scan direction and a scan azimuth of 48.63°

- Image 2 has an incidence of -28.34° in the scan direction, 15.28° in the ortho-scan direction and a scan azimuth of 24.34°

`howstereo.py --inc1 -5.44 7.07 --az1 48.63 --inc2 -28.34 15.28 --az2 24.34`

Output:

```
pair	b/h	angle
im1-im2	0.48	27.1°
```

#### Example 2:

The information for each image is listed in image_list.txt:
```
15,10,120
-15,-10,120
20,0,30
-20,0,30
```

`howstereo.py --input_file image_list.txt --show_plot`

Output:

```
pair	b/h	angle
im1-im3	0.32	18.0°
im2-im4	0.32	18.0°
im1-im4	0.60	33.2°
im2-im3	0.60	33.2°
im1-im2	0.64	35.7°
im3-im4	0.73	40.0°
```

<img src="https://github.com/IPGP/image-pairs-stereo-capacity/blob/master/Figure.jpg" alt="Interactive plot of the result (screen capture)" width=400>

#### Example 3:

In this example we are comparing, for two SPOT6 images, the B/H found using 1) the target azimuth from the DIMAP files and 2) the scan azimuth from the GeoStore.

- Image 1 has an incidence of -19.42° in the scan direction, -0.27° in the ortho-scan direction and a **target** azimuth of 2.11°

- Image 2 has an incidence of 17.16° in the scan direction, -9.14° in the ortho-scan direction and a **target** azimuth of 203.79°

`howstereo.py --inc1 -19.42 -0.27 --az1 2.11 --inc2 17.16 -9.14 --az2 203.79 --az_mode target`

Output:

```
pair	b/h	angle
im1-im2	0.69	38.0°
```

Which is equivalent to:

- Image 1 has an incidence of -19.42° in the scan direction, -0.27° in the ortho-scan direction and a **scan** azimuth of 182.88°

- Image 2 has an incidence of 17.16° in the scan direction, -9.14° in the ortho-scan direction and a **scan** azimuth of 176.27°

`howstereo.py --inc1 -19.42 -0.27 --az1 182.88 --inc2 17.16 -9.14 --az2 176.27 --az_mode scan`

Output:

```
pair	b/h	angle
im1-im2	0.69	38.0°
```
