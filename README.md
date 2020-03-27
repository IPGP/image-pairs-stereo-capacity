# image-pairs-stereo-capability
Estimate the stereoscopic capability (B/H) of image pairs from the Pleiades and SPOT6-7 satellites

### Usage:

```
howstereo.py --incid1 scan ortho --az1 azimuth --incid2 scan ortho --az2 azimuth [-h] [--show_plot]
```
Mandatory arguments:
```
--incid1 scan ortho  incidence angles for image 1 (in degrees)
--az1 azimuth        scan azimuth for image 1 (in degrees)
--incid2 scan ortho  incidence angles for image 2 (in degrees)
--az2 azimuth        scan azimuth for image 2 (in degrees)
```
Optional arguments:
```
-h, --help           show this help message and exit
--show_plot          show a 3D, interactive plot
```

### Example:
- Image 1 has an incidence of -5.44° in the scan direction, 7.07° in the ortho-scan direction and a scan azimuth of 48.63°.
- Image 2 has an incidence of -28.34° in the scan direction, 15.28° in the ortho-scan direction and a scan azimuth of 24.34°

`howstereo.py --incid1 -5.44 7.07 --az1 48.63 --incid2 -28.34 15.28 --az2 24.34 --show_plot`
