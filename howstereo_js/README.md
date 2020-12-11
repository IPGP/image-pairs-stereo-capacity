# javascript version

## Use

### Insert in your page header

``` html
 <script src="url/to/howstereo.js"></script>
```

### Or from cdn

``` html
<script src="https://cdn.jsdelivr.net/gh/IPGP/image-pairs-stereo-capacity@[tag]/howstereo_js/howstereo.js"></script>
```

### Use in your script

``` js
  // create 2 images from scan, ortho and azimuth data in degree
  var img1 = howstereo.img(-5.44, 7.07, 48.63)
  var img2 = howstereo.img(-28.34, 15.28, 24.34)
  // compute b_to_h 
  var hs = howstereo.compute_b_to_h(img1, img2)
  // check results
  console.log(hs.angle)
  // display 27.1
  console.log(hs.bh)
  // display 0.48

```

