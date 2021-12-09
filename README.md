# `ece4522-digital_signals_processing`
MATLAB projects for Digital Signals Processing (ECE 4522).

## [`/assignment2/`](assignment2)

Introduces various types of systems and provides hands on experience
with the impulse response. It cements the concept of convolution using
MATLAB including effects on a sound sample, and expands it with the
concept of a restoring system to undo convolution. We will find the
error plot and worst-case error to compare out actual result with the
expected input and prove that a restoring system is possible.

### Terms

convolution, finite-duration impulse response, FIR, discrete-time
system, echo system, restoration system, digital signal processing,
cascading systems

## [`/assignment3/`](assignment3)

Introduces the backward difference system, and uses cascading to
produce an edge detector. This works because the difference of the
current value and the previous are zero if they were equal, and larger
if they are more different. The backward difference system is produced
in 3 steps, and is first tested. We finally apply the edge detector to
a photograph of a barcode to prepare it for decoding.

### Terms

convolution, filter, backward difference system, difference system,
edge detection, threshold, barcode, encode, decode

## [`/doc/`](doc)

The source code for the reports of each assignment in
L<sup>A</sup>T<sub>E</sub>X. The IEEE Conference Template (2019)<!--
--><sup>[&lsqb;1&rsqb;](#ref-1)</sup> is used.

## References

1. <a name='ref-1'></a> IEEE. (2019). *IEEE Conference Template*. Overleaf. Retrieved from
<<https://www.overleaf.com/latex/templates/ieee-conference-template/grfzhhncsfqn>>

