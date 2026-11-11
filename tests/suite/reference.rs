//! Results of the Python implementation (version 0.1.0.dev0) for the public
//! reference machines, used as regression values. Generated from its output
//! files; do not edit by hand.

/// Configuration values of the `w7x` case.
pub(crate) const W7X: [(&str, f64); 34] = [
    ("symmetry", 5.0),
    ("coilspermodule", 10.0),
    ("rmajor_ref", 5.511_939_539_133_372_5),
    ("rminor_ref", 0.512_692_195_535_189_6),
    ("aspect_ref", 10.750_972_195_665_206),
    ("bt_ref", 3.000_000_000_000_000_4),
    ("vol_plasma", 28.598_786_068_665_905),
    ("plasma_surface", 136.662_192_597_402_22),
    ("coil_rmajor", 5.617_761_381_926_321),
    ("coil_rminor", 1.144_289_593_869_747_3),
    ("coillength", 425.862_628_660_959_73),
    ("coilsurface", 300.637_175_675_506_5),
    ("maximal_coil_height", 3.096_473_429_003_781),
    ("dmin", 0.279_206_983_017_383),
    ("max_portsize_width", 0.768_094_341_319_239),
    ("min_bend_radius", 0.395_593_474_001_869_8),
    ("min_plasma_coil_distance", 0.393_692_865_778_217_8),
    ("derivative_min_LCFS_coils_dist", -0.435_858_665_116_085_3),
    ("WP_area", 0.027),
    ("WP_ratio", 1.2),
    ("i0", 1.745_239_494_739_875),
    ("WP_bmax", 6.292_455_758_602_379_5),
    ("a1", 0.640_558_303_072_252),
    ("a2", 0.031_125_225_355_697_214),
    ("inductance", 0.000_418_194_629_645_164_8),
    ("max_force_density", 153.632_917_914_563_62),
    ("max_force_density_MNm", 4.148_088_783_693_218),
    ("max_lateral_force_density", 128.239_657_427_817_4),
    ("max_radial_force_density", 147.061_189_942_686_87),
    ("centering_force_max_MN", 2.910_502_480_053_418),
    ("centering_force_min_MN", -1.286_227_062_388_286_5),
    ("centering_force_avg_MN", 1.551_325_570_978_960_8),
    ("epseff", 0.009_5),
    ("neutron_peakfactor", 1.6),
];
/// Surface lofted between the coils of the `w7x` case [m^2].
pub(crate) const W7X_LOFTED: f64 = 292.425_980_478_610_9;
/// Magnetic energy of the `w7x` case [GJ].
pub(crate) const W7X_ENERGY: f64 = 0.636_881_334_258_488_3;

/// Configuration values of the `hsx` case.
pub(crate) const HSX: [(&str, f64); 34] = [
    ("symmetry", 4.0),
    ("coilspermodule", 12.0),
    ("rmajor_ref", 1.211_301_676_144_077),
    ("rminor_ref", 0.121_493_658_017_303_5),
    ("aspect_ref", 9.970_081_532_746_011),
    ("bt_ref", 1.003_439_605_182_259_6),
    ("vol_plasma", 0.352_930_568_622_345_06),
    ("plasma_surface", 7.758_471_556_742_802),
    ("coil_rmajor", 1.237_427_231_458_200_4),
    ("coil_rminor", 0.316_838_430_417_726_14),
    ("coillength", 107.396_730_267_511_09),
    ("coilsurface", 17.395_998_413_667_453),
    ("maximal_coil_height", 1.134_473_535_940_21),
    ("dmin", 0.092_901_056_241_991_3),
    ("max_portsize_width", 0.193_647_508_895_788_2),
    ("min_bend_radius", 0.082_210_867_435_301_61),
    ("min_plasma_coil_distance", 0.142_915_829_267_282_05),
    ("derivative_min_LCFS_coils_dist", -0.539_668_349_812_674_5),
    ("WP_area", 0.003),
    ("WP_ratio", 1.2),
    ("i0", 0.150_072_550_000_000_03),
    ("WP_bmax", 1.905_705_587_979_287_2),
    ("a1", 0.601_683_145_876_238_6),
    ("a2", 0.029_651_103_952_555_03),
    ("inductance", 0.000_101_726_938_960_008_13),
    ("max_force_density", 40.086_401_753_617_67),
    ("max_force_density_MNm", 0.120_259_205_260_853_01),
    ("max_lateral_force_density", 37.305_247_515_572_53),
    ("max_radial_force_density", 38.850_286_362_829_98),
    ("centering_force_max_MN", 0.040_498_445_429_084_79),
    ("centering_force_min_MN", -0.013_659_272_769_521_62),
    ("centering_force_avg_MN", 0.012_176_972_410_738_379),
    ("epseff", 0.005),
    ("neutron_peakfactor", 1.6),
];
/// Surface lofted between the coils of the `hsx` case [m^2].
pub(crate) const HSX_LOFTED: f64 = 18.708_103_717_934_06;
/// Magnetic energy of the `hsx` case [GJ].
pub(crate) const HSX_ENERGY: f64 = 0.001_145_535_374_433_323;

/// Configuration values of the `w7x_scaled` case.
pub(crate) const W7X_SCALED: [(&str, f64); 34] = [
    ("symmetry", 5.0),
    ("coilspermodule", 10.0),
    ("rmajor_ref", 19.999_999_999_999_97),
    ("rminor_ref", 1.860_296_876_971_184),
    ("aspect_ref", 10.750_972_195_665_183),
    ("bt_ref", 5.500_000_000_000_002),
    ("vol_plasma", 1_366.231_362_904_314),
    ("plasma_surface", 1_799.283_035_837_789),
    ("coil_rmajor", 20.383_973_162_410_95),
    ("coil_rminor", 4.152_039_715_768_22),
    ("coillength", 1_545.236_937_515_163_2),
    ("coilsurface", 3_958.163_994_402_402_7),
    ("maximal_coil_height", 11.235_513_042_258_921),
    ("dmin", 1.013_098_859_430_819_5),
    ("max_portsize_width", 2.787_020_198_120_710_6),
    ("min_bend_radius", 1.435_405_708_619_479),
    ("min_plasma_coil_distance", 1.428_509_376_719_725_2),
    ("derivative_min_LCFS_coils_dist", -0.435_858_665_116_294),
    ("WP_area", 0.3),
    ("WP_ratio", 1.2),
    ("i0", 11.609_727_275_272_338),
    ("WP_bmax", 12.145_617_878_790_338),
    ("a1", 0.638_398_882_474_363_6),
    ("a2", 0.031_171_753_278_231_1),
    ("inductance", 0.001_543_361_191_717_309_6),
    ("max_force_density", 169.675_601_907_230_7),
    ("max_force_density_MNm", 50.902_680_572_169_22),
    ("max_lateral_force_density", 142.572_087_754_704_9),
    ("max_radial_force_density", 160.123_344_961_677_45),
    ("centering_force_max_MN", 128.958_375_810_814_06),
    ("centering_force_min_MN", -57.013_257_558_380_2),
    ("centering_force_avg_MN", 68.672_001_796_582_85),
    ("epseff", 0.009_5),
    ("neutron_peakfactor", 1.6),
];
/// Surface lofted between the coils of the `w7x_scaled` case [m^2].
pub(crate) const W7X_SCALED_LOFTED: f64 = 3_850.056_082_909_638_5;
/// Magnetic energy of the `w7x_scaled` case [GJ].
pub(crate) const W7X_SCALED_ENERGY: f64 = 104.011_561_305_284_37;
