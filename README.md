# High precision subsample shift detection and correction in Odin
This is what appears to be a very precise shift detection.

## Description
I was searching for ways to increase the accuracy of a shift detection by cross correlation implemented with FFT and with with phase detection, and after many iterative steps with a llm and after many manual correction and modificationn to the code. I come up with this very precise implementation ( low error ) that I share in this repository.

## License
MIT Open Source license

## Example of output and precision

```bash
(base) joaocarvalho@soundofsilence:~/zed_editor/shift_1d_deteection> time ./shift_1d_detection.exe

Begin subsample vector 1D shift detection and correction ...

N : 8192

Integer shift from cross-correlation peak:  -15 samples
Calculated subsample shift :                -15.1234567889414393 samples
True shift :                                -15.1234567890123461 samples
Error :                                     7.0906835958339798e-11 samples
Calculated subsample shift ( after shift ): -0.0000000001314220 samples

Integer shift from cross-correlation peak:  19 samples
Calculated subsample shift :                19.1234567888950551 samples
True shift :                                19.1234567890123444 samples
Error :                                     1.1728928939191974e-10 samples
Calculated subsample shift ( after shift ): 0.0000000000209184 samples

Integer shift from cross-correlation peak:  100 samples
Calculated subsample shift :                100.1234567889132450 samples
True shift :                                100.1234567890123515 samples
Error :                                     9.9106500783818774e-11 samples
Calculated subsample shift ( after shift ): 0.0000000001664375 samples

Integer shift from cross-correlation peak:  -100 samples
Calculated subsample shift :                -100.1234567890396647 samples
True shift :                                -100.1234567890123515 samples
Error :                                     2.7313262762618251e-11 samples
Calculated subsample shift ( after shift ): 0.0000000000618456 samples

Integer shift from cross-correlation peak:  4000 samples
Calculated subsample shift :                4000.1234567889960090 samples
True shift :                                4000.1234567890123799 samples
Error :                                     1.6370904631912708e-11 samples
Calculated subsample shift ( after shift ): 0.0000000000800355 samples

Integer shift from cross-correlation peak:  -4000 samples
Calculated subsample shift :                -4000.1234567889787286 samples
True shift :                                -4000.1234567890123799 samples
Error :                                     3.3651303965598345e-11 samples
Calculated subsample shift ( after shift ): -0.0000000000927685 samples


 ... end subsample vector 1D shift detection and correction.


real	0m0,122s
user	0m0,112s
sys	0m0,010s
```

## Have fun
Best regards, <br>
Joao Carvalho
