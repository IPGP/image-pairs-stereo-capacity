Number.prototype.toRad = function () {
  return this * Math.PI / 180
}
Number.prototype.toDeg = function () {
  return this * 180 / Math.PI
}
const howstereo = {
  img (scan, ortho, az) {
    return {
      scan: scan,
      ortho: ortho,
      az: az,
      s_comp: Math.cos(ortho.toRad()) * Math.sin(scan.toRad()),
      o_comp: Math.cos(scan.toRad()) * Math.sin(ortho.toRad()),
      z_comp: Math.cos(ortho.toRad()) * Math.cos(scan.toRad())
    }
  },
  compute_b_to_h (im1, im2) {
    var tetaDeg = im2.az - im1.az
    var teta = tetaDeg.toRad()
    var p = [im1.s_comp, im1.o_comp, im1.z_comp]
    var p_prime = [
      Math.cos(teta) * p[0] - Math.sin(teta) * p[1],
      Math.sin(teta) * p[0] + Math.cos(teta) * p[1],
      p[2]
    ]
    var q = [im2.s_comp, im2.o_comp, im2.z_comp]
    var stereo_angle = Math.acos(p_prime[0] * q[0] + p_prime[1] * q[1] + p_prime[2] * q[2])
    var half_angle = stereo_angle / 2.0
    var b_to_h = Math.tan(half_angle) * 2
    return {angle: stereo_angle.toDeg(), bh: b_to_h}
  }
}