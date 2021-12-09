# `ece4522-digital_signals_processing`
MATLAB projects for Digital Signals Processing (ECE 4522).

## <a name='assignment2'></a> [`/assignment2/`](assignment2)

Introduces various types of systems and provides hands on experience
with the impulse response. It cements the concept of convolution using
MATLAB including effects on a sound sample, and expands it with the
concept of a restoring system to undo convolution. We will find the
error plot and worst-case error to compare out actual result with the
expected input and prove that a restoring system is possible.

### <a name='assignment2-terms'></a> Terms

convolution, finite-duration impulse response, FIR, discrete-time
system, echo system, restoration system, digital signal processing,
cascading systems

## <a name='assignment3'></a> [`/assignment3/`](assignment3)

Introduces the backward difference system, and uses cascading to
produce an edge detector. This works because the difference of the
current value and the previous are zero if they were equal, and larger
if they are more different. The backward difference system is produced
in 3 steps, and is first tested. We finally apply the edge detector to
a photograph of a barcode to prepare it for decoding.

### <a name='assignment3-terms'></a> Terms

convolution, filter, backward difference system, difference system,
edge detection, threshold, barcode, encode, decode

## <a name='doc'></a> [`/doc/`](doc)

The source code for the reports of each assignment in
L<sup>A</sup>T<sub>E</sub>X. The IEEE Conference Template (2019)<!--
--><sup>[[1]](#references)</sup> is used.

## <a name='references'></a> References

1. IEEE. (2019). *IEEE Conference Template*. Overleaf. Retrieved from
<https://www.overleaf.com/latex/templates/ieee-conference-template/grfzhhncsfqn>

