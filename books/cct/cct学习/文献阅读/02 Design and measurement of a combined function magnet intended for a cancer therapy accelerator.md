# 02 Design and measurement of a combined function magnet intended for a cancer therapy accelerator

https://www.researchgate.net/publication/32154504_Design_and_measurement_of_a_combined_function_magnet_intended_for_a_cancer_therapy_accelerator

时间 2001 作者 A Morita，Y Iwashita，A Noda，T Shirai，M Tadokoro 机构 日本日立 Hitachi

### abstract

A compact proton synchrotron using combined function magnets is proposed to help realize the wider
availability of charged particle cancer therapy facilities. This combined function magnet was designed
with the help of three-dimensional magnetic field calculations to take account of a realistic fringe and the
interference among the magnetic poles. An evaluation scheme for tune values based on particle tracking
was developed to improve the magnet design. To verify the magnet design, a model magnet was fabricated
and measured. In order to achieve a tune value evaluation from the measured magnetic field, schemes
for accurate field mapping and field interpolation were developed. From the tune value evaluation of the
measured magnetic field, it was thought that the performance of the model magnet was good enough to
construct a synchrotron. In this paper, we report details of the design and the evaluation scheme for the
combined function magnet and the results of the field measurements of the model magnet.

### I. INTRODUCTION

Radiation cancer therapy has recently been the subject
of attention because of both its ability to preserve human
body function and shape and its light load to the patient
compared with other therapies. Among the radiation therapies, charged particle therapy has the advantage that it
can localize dose distribution to the tumor, largely due to
the presence of the Bragg peak which reduces damage to
normal cells. But such charged particle cancer therapy is
not yet common because of the high cost of the accelerator.


To realize widespread use of charged particle cancer
therapy, an accelerator with low construction cost and easy
handling is required. A compact proton synchrotron using
combined function magnets is proposed as one approach
for such an accelerator system [1,2].

Using combined function bending magnets as the main
focusing elements, tracking between bending dipole and
focusing quadrupole becomes unnecessary because the
functions of bending and focusing are both realized with
the same magnets. The operating point of the ring will
be fixed if the relative distribution of the magnetic field
does not change over the whole excitation range. As a
merit of the combined function magnet, the possibility of
eliminating a tune change due to a power supply ripple
is pointed out and an experiment of proof of principle
was achieved [3,4]. Power, coil space, and alignment
costs between dipole and quadrupole can be saved. Using
an untuned cavity [5,6] and an rf source controlled by
magnetic field strength, the acceleration operation will
be easy.

However, the adjustment of the tune values after
fabrication of such a magnet is not easy. To achieve a
good magnet design we employed the three-dimensional
magnetic field calculation code. To evaluate and adjust
the tune values of the magnet, an evaluation method by
particle tracking using calculated field distribution was
developed.

In order to verify this design scheme, a model magnet
was fabricated and its magnetic field distribution was measured. This paper describes our design of a combined function synchrotron emphasizing the design of the magnet and
the measured results of the model magnet.

### II. DESIGN

### A. Initial design

The synchrotron ring is composed of six sets of 2 m
drift space and 60± sector magnets (Fig. 1). In the orbit in
the bending magnet, the radius of curvature and the maximum field strength are designed as 1.9 m and 1.28 T, respectively, for protons up to 240 MeV. The bending magnet has three combined function magnet sections and an
FDF triplet focusing structure. The bending angles of the
three poles are designed to be 15±, 30±, and 15±. The design target values for horizontal and vertical tunes are both
1.75. After a preliminary calculation, the n indices of the
magnet poles [n = -(ρ/B)(∂B/∂ρ)] were determined

<img src="./img/002-01.jpg"></img>

as 25.855 for the F sector and 6.164 for the D sector. The
main coil is made by a hollow conductor, and the winding
number of the main coil is 28 per pole. The conductors
of the main coil make a 7 * 4 * W * H conductor stack,
and this is installed into the magnetic pole with a one turn
correction coil.

The basic magnet pole shapes are designed to suppress
the sextupole component by a two-dimensional magnetic
potential formula,

<img src="./img/002-02.jpg"></img>

where r, g, and n are the radius of the curvature of the
magnet sector, the gap height on the designed orbit, and
the n index of the magnet pole, respectively. x and y
denote the horizontal and vertical distance, respectively,
from the designed orbit. The gap height on the designed
orbit is 67 mm. The sextupole error field is evaluated from
a magnetic field distribution calculated by 2D code and
its polynomial fitting. In order to avoid a variation of
the field distribution due to saturation, slots to equalize
the magnetic flux density [7] are situated at the low flux
density area in the pole, as shown in Fig. 2. From the
magnetic field calculations of the initial design, a useful
aperture, +/-70 mm, of the magnet with an error of 0.1%
was obtained.

The entrances of the sectors are shaped by a step function which approximates the magnetic potential of the
Rogowski cut pole. There are five steps. The width of
these Rogowski-like steps are 10, 10, 10, 20, and 20 mm,
sequentially, from the edge of the sectors. Figure 3 shows
the original Rogowski curve and the step function of the
real pole cut.

<img src="./img/002-03.jpg"></img>

<img src="./img/002-04.jpg"></img>

### B. Operating point evaluation method

To evaluate a realistic magnetic field including fringe
fields and interference between the F and D sectors, the
magnetic field distribution of the model magnet was calculated by the 3D magnetic field calculation code, TOSCA.
A conventional analysis which used the transfer matrices
of the ring component was not enough to evaluate the tune
values from such a calculated magnetic field. Therefore
we introduced a tune value evaluation method based on
particle tracking [8]. In this method, the condition of the
reference particle which made the closed orbit is determined from the given magnetic field at the first step. At
the second step, a curvilinear coordinate s, x, y, which is
based on the trajectory of the reference particle, is introduced and becomes the basis of the trajectory description.
At the third step, a set of particles that has an initial distribution in a phase space is tracked in the magnet. Finally,
the transfer matrices are reconstructed from the initial and
final phase space distributions using least squares fitting.


From the curvilinear coordinate, the development of the
trajectory is considered as the mapping of the phase space,

m(s1|s0):X0 -> X1

where s0 and s1 are the curvilinear coordinates of the initial
and the mapped points, respectively, and x0 and x1 are
the initial phase space vector at the point of s0 and the
mapped vector at the point of s1, respectively. The linear
component is the lowest order of this mapping because
the origin is the fixed point of mapping from the definition

of the reference orbit of the coordinate. In this mapping,
the transfer matrix is understood as a linear representation.
To extract the transfer matrix Ms1js0, we reconstructed
it from the two phase space distributions by the following
method. The initial vector of the particle number n is x
n is X0 and the mapped vector of n is x
n is X1 ; the least squares method
for the transfer matrix is then written as

<img src="./img/002-05.jpg"></img>

where T denotes the transposition of the vector and the matrix. This reconstructed transfer matrix Ms1js0 is equal to the linear representation of the mapping ms1js0 under the L2 norm.


Once a transfer matrix of an arbitrary interval of the ring is calculated, we can obtain a transfer matrix of one revolution.
The phase advance of betatron oscillation m and the Twiss parameter a, b, and g are easily obtained from the twodimensional transfer matrix of one revolution using the periodic boundary condition of the betatron oscillation in the following relationship:

<img src="./img/002-06.jpg"></img>

where Mrevs, as, bs, and gs are the transfer matrix of one revolution begun from the point s and Twiss
parameters at the point s, respectively. mrev denotes the
phase advance of the betatron oscillation per revolution.
Finally, the tune value n is derived from the phase advance
per revolution mrev as

v = μ(rev)/2π

When we applied this method to our magnet, it was assumed that the closed reference orbit stayed on the median
plane and coincided with the designed orbit at the middle
of the drift space. With the help of these assumptions, we
define the initial position and the initial direction of the
reference orbit. The last parameter needed to complete the
definition of the closed reference orbit is the initial momentum of the orbit. Therefore, a candidate of the closed
reference orbit corresponding to the given initial momentum can be obtained by tracking the particle trajectory.
Under good field symmetry, the closed reference orbit parallels the designed orbit at the center of the magnet. Thus,
assuming such good field symmetry, the closed reference
orbit is obtained by choosing the trajectory that parallels
the designed orbit at the center of the magnet from the candidates generated by the tracking method. The angle between the designed orbit and the tracked trajectory at the
center of the magnet is a function of the initial momentum of the trajectory. By applying the Newton-Raphson
method to this function, the initial momentum generating
the closed reference orbit is obtained. Once the initial
momentum of the closed reference orbit is obtained, the
curvilinear coordinate with respect to the reference orbit is
defined and the calculations on the curvilinear coordinate
can be performed. In order to make it easy to calculate the
inverse matrices, the initial distribution of xn0 was selected

to eliminate the off diagonal terms of [xxx]
Eq. (3). In the evaluation of the Twiss parameters and
the phase advance, it was assumed that the betatron oscillations could be decoupled into x-x0 and y-y0 subphase
spaces. This treatment, to decouple the betatron oscillations, is valid under the assumption of the magnet having
good symmetric poles.

### C. Evaluation and tuning of the rough design

The operating point in our rough design is evaluated as
nh, ny1.64, 1.86. These tune values can be adjusted
by modifying the ratio of bending angles between the F
and D sectors. This increases the bending angle of the F
sector by Du degrees and decreases the bending angle of
the D sector by 2Du degrees, as shown in Fig. 4. Because
we plan to make the magnetic pole of the laminated iron
sheets, this modification can be easily incorporated into the
design by changing the number of the iron sheets laminated
of the F and D sectors. Changing the ratio of the bending

<img src="./img/002-07.jpg"></img>

<img src="./img/002-08.jpg"></img>

angles can be treated as a local modulation of the n index
around the transition region between the F and D sectors.
From the tune shift formula of Hill’s equation, the tune
shift per bending angle modulation is obtained as follows:

<img src="./img/002-09.jpg"></img>

where nF and nD are design n indices and r is the radius
of the designed orbit curvature. bh and by are horizontal and vertical beta functions at the border between the
F and D sectors, respectively. Thus, under the linear approximation, the ratio of the horizontal and vertical tune
shifts is given by the ratio of the beta function around the
transition region, and the magnitude of the tune shifts is
proportional to the modulation of the bending angle Du.
Figure 5 shows the beta functions evaluated from the calculated magnetic field distribution, where bh and by in
Eqs. (6) and (7) are about 2 and 4 m, respectively. Thus the
estimated proportional coefficients of the tune shifts for the
angle Du are 0.21 deg21 horizontally and 20.42 deg21

vertically. The tune shifts predicted by Eqs. (6) and (7)
are shown in Fig. 6 by the dotted line from the open circle
labeled “0.0 deg.” The open and the filled circles in Fig. 6
show the evaluated tune values by the tracking method.
These circles have good agreement with the predicted line.
Therefore, using this Du tuning method, we can control the
operating point on the predicted line, as shown in Fig. 6.

### D. Final design and model magnet

Considering the operating point and resonance lines in
Fig. 6, we selected 0.25± for the Du parameter and operating point 1.70, 1.74 of the final point to avoid resonance
lines up to the fifth order. Figure 7 shows the excitation dependence of tune values evaluated from the magnetic field
calculated by TOSCA with main coil current conditions of
383, 574, 765, 861, 957, 1053, 1148, 1244, and 1349 A.
In Fig. 7, the horizontal tune is shifted to the lower side
with an excitation increase, and the vertical tune has a peak
around current I  1000 A. From the analysis of the position dependence of the contributions to the tune shift,
the qualitative property of excitation dependence can be
explained as follows. The ratio between the vertical and
horizontal beta functions, which determine the modulation
magnitude of the tune values, is especially large in the D
sector. In the low excitation case, field deformation in the
D sector is small, and the tune shifts are dominated by field
deformation of the F sector. Then the vertical and horizontal tunes move in opposite directions. In the condition beyond the current of 1000 A, a strong field deformation due
to pole saturation appears in the D sector. This saturation
affect dominates only the vertical tune, then the direction
of the vertical tune shift is changed. Therefore it is easily
supposed that the tune excursion curve is modified by the
B-H curve of iron, the lamination packing factor, and the
machining error of the slots equalizing flux density, etc.

<img src="./img/002-10.jpg"></img>

To verify the three-dimensional calculation, we fabricated a model magnet reflecting the final design and measured the field distribution. Figure 8 shows the lower
half of the model magnet being assembled. The poles
of the magnet are made of laminated silica steel sheets
of 0.5 mm thickness compatible with Nippon Steel Corporation 50H600 (see http://www.nsc.co.jp/si-steel/products/
05.html). The achieved packing factor of the laminated
pole was about 0.95.




<img src="./img/002-11.jpg"></img>

### III. FIELD MEASUREMENT
### A. Basic strategy

In order to check the n index and the operating point of
the model magnet, we needed accurate field gradient additions to the field strength. Because a harmonic coil method
is not feasible for a sectored magnet, Hall probe mapping
was chosen. This Hall probe mapping was achieved using the combination of a three-axes stage and Hall probes.
The alignment method and the details of the measurement
setup are described in a later section.


The information of field distribution required by the
transfer matrix reconstruction using particle tracking is less
than that in a full three-dimensional distribution, because
of its linear approximation, magnet pole symmetry, and
Maxwell law. Although evaluation of the dynamic aperture requires a full three-dimensional measurement, such a
measurement is difficult because of the limitations of the
measurable area, the probe stability, and the measurement
time. Thus, we restricted measurement to the minimum
set, which was enough to evaluate the n index and tune
values by tracking. The symmetry between the upper and
lower poles and Maxwell equation ∇ B = 0 told us that
the magnetic field flux crosses perpendicular to the median
plane, and that the magnetic field near the median plane
can be reconstructed from the major component Bz on the
median plane as follows:


<img src="./img/002-12.jpg"></img>

### B. Alignment method

Before field mapping, both the three-axes stage and the
Hall sensor axis of the Hall probe have to be mechanically
aligned against the magnet and the magnetic median plane
found. To arrange the horizontal plane of the stage with
the horizontal plane of the magnet, both the magnet and
the stage are leveled using a water level. To arrange the
movable area and the moving axes of the stage within
the measurement region, the stage axis is aligned to the
magnet axis at an appointed angle using a theodolite in
the horizontal plane after translation of the stage. After
mechanical stage alignment, the probe axes have to be
aligned with the coordinate axes of the magnetic field to
accurately measure the component of the magnetic field. It
is difficult to find the coordinate axes of the magnetic field
because of the change in the direction of the magnetic flux
caused by the combined function magnet. Thus, a plane
in which the direction of the flux is already known has to
be found before alignment of the probe axes. Fortunately,
the median plane has such a property because the whole
magnetic flux is perpendicular to the median plane if the
poles are symmetric. Because the flux direction in the offmedian plane depends on the n index of the pole, the
positional dependence of the flux direction over the whole
measurement region can be detected by a roughly aligned
probe. Thus the problem of finding the magnetic median
plane is replaced with one of finding the plane where the
directions of the flux of the F and D sectors are aligned.
Using this searching method for the median plane, the
accuracy requirement of the primary probe alignment can
be relaxed. Once the median plane is found, the probe axes
can be aligned with that plane and the alignment error of
the probe size effect reduced by iterating the median plane
searching and probe alignment. Finally, the probe origin in
the horizontal plane is determined mechanically by using
the geometry of the probe and the magnet.

### C. Setup and equipment

It is difficult to measure a field map of a whole area with
a single setup because of the sector shape of the magnet.
Considering the mirror symmetry of the magnet along the
beam path, we measured only half of the magnet area with
a small overlap of the mirrored area. Figure 9 shows the
geometry of the measurement setup. In this setup, the
probe arm has a tilt angle of 15± to the central orbit at
the entrance, to cover the maximum magnet aperture of





<img src="./img/002-13.jpg"></img>

the measurement area. The longest axis of the movable
stage is aligned to the probe arm direction, which is used
as the main scan axis to reduce the vibration of the probe
arm when excited by the acceleration and deceleration of
the stage. Figure 10 shows the setup as seen from the rear
of the stage.
Measurement was performed on 5 * 5 mm rectangular
grids along the stage axes. Considering the requirements
for the evaluation of the field characteristics, only the field
distribution around the designed orbit is needed. The zonal
area of the field distribution around the closed orbit, whose
width is proportional to the beta function, is required to
track the betatron particles. The beta function is maximized at the center of the D sector, and this position stays
at the deepest area in the measurement region. Hence, the

scan width of this measurement area cannot be extended
widely because of the conflict between the probe arm and
the coil support. In order to reduce the measurement time,
we decided to cut the measurement points which had relatively small contributions to the tracking aperture and the
characteristic evaluation. Thus the selected measurement
area covered half of the D sector, the F sector, and a
10± fringe within a radius from 1.84 to 1.96 m. The total number of the measurement points in the selected area
was 7139.


Our method for searching the median plane requires
a two- or three-axes Hall probe. Two commercially
available three-axes Hall probes could not satisfy our
requirements for the measurement stability, reproducibility, and differential linearity. We finally decided to use
a single Hall sensor Group 3 MPT-141 that has a small
sensitive area, high enough accuracy, and good stability (see http://ourworld.compuserve.com/homepages/
group3tech/DTM.htm). Its sensitive area is a rectangle of
1.0 3 0.5 mm. The maximum temperature coefficient of
the corrected readout is guaranteed within 610 ppm±C
according to the catalog specification. Two Hall MPT-141
sensors were mounted on the aluminum mount (see
Fig. 11), in which the geometric center of the sensing
areas of the two Hall sensors stayed on the same horizontal
plane at a distance of 7.5 mm. This probe mount was
attached to the tip of the probe arm. To control the Hall
sensor angle against the flux, the directions of the two Hall
sensors were perpendicular to the probe arm. The probe
arm was fixed on the three-axes stage, which was driven
with a resolution of 5 mm by three stepping motors. The
position reproducibility of the motor system, measured
by a linear scale with 5 mm accuracy, was held within
20 mm on the axis having the worst reproducibility. The
stepping motor drivers of the three-axis stage were connected to a TUJI DENSHI PM4C four-channel pulse motor controller (see http://www.tsuji-denshi.co.jp/english/
pm4c05a.html). All measurement devices were controlled
by a Note PC via a NI PCMCIA-GPIB card (see http://
www.ni.com/catalog/pdf/1gpib756a.pdf). Even though
the linear scalers installed in each axis were independent
of the control system, it was useful to confirm that the
motor system worked well.


<img src="./img/002-14.jpg"></img>

### D. Field interpolation

Evaluation of the field characteristics and particle tracking required a continuous smooth field distribution. The
Fourier series expansion was chosen for interpolation because of its easy handling of derivatives. Using a Cartesian
measurement grid, a Fourier series of a field distribution
was obtained by two-dimensional discrete cosine transformation (2D DCT). Although all the measured points were
on rectangular grid points, as previously described, not all
grid points in the rectangular area had data because of the
sector shape of the magnet. A grid point that had no measured data was treated as a free parameter with an initial
value of zero. The following iterations were applied to reduce the number of Fourier components and noise. First,
the Fourier components were obtained by DCT. At the
second step, higher frequency Fourier components were
reduced by a low pass filter (LPF). At the third step, the
field map was reconstructed by inverse DCT (IDCT), including the missing data points. At the final step, the interpolation error was estimated on the measured grid points,
and the calculated data in the measured area were replaced
by measured values. The sequence DCT, LPF, IDCT, and
refilling were iterated until the interpolation error came
within tolerances.

### IV. RESULTS OF MEASUREMENTS

### A. Characteristics of the model magnet

The excitation curve was initially compared with the
TOSCA calculations. The excitation curve at the center of
the D sector is shown in Fig. 12. The open and the closed
circles show the results of magnetic field calculations and
the results of field measurements, respectively. The excitation curve of the measurements almost corresponds to
the curve of the field calculations. The apparent threshold
of the saturation phenomenon exists around the main coil
current of I = 1200 A. But the nonlinearity caused by
magnetic pole saturation can be seen above a coil current
of I = 900 A from the curve of the field strength normalized by the coil current.




<img src="./img/002-15.jpg"></img>


<img src="./img/002-16.jpg"></img>


<img src="./img/002-17.jpg"></img>


In order to compare characteristics of field distribution,
the field distributions of TOSCA calculations were compared
with the measured ones by evaluating the magnetic field
strength normalized by the main coil current BI and the
n index XXX distribution in the designed orbit.


Figures 13 and 14 show the normalized field strengths
and n indices, respectively, on an arc of r  1.9 m. This
arc overlaps with the designed orbit of the magnet. Although the position of the magnetic field edge and the field
strength for low field cases are well predicted by the calculation, the field strength of the fringe and the absolute
value of the n index are different.


The lamination structure of the pole edges was different between the model magnet and the TOSCA calculation
model. In the TOSCA calculation, the whole magnet is constructed by sector lamination. But the Rogowski-like pole
edge of the model magnet is laminated linearly because of
assembling restrictions. This difference of the edge structure may be one of the causes of the characteristic difference around the magnet edges.

### B. Tune diagram of the measured result

Previous graphs have shown the fidelity of the model
magnet against the proposed design. The information is
not enough to verify the operating point of the whole ring
because of closed orbit distortion. Thus, we have to evaluate tune values by a tracking method based on the measured
magnetic field distribution. Figure 15 shows an operating
point evaluated by the tracking method. The field of the
same excitation current was measured two times with a few
days separation. The open circles show the primary measurement results and the filled circles show the secondary
measurement results. These open and filled circles show
the reproducibility of the measurements. The tune difference between the two measurements is larger in the lower
field case because of relative error emphasis.


The major operating point was evaluated from the measurement as nh, ny1.71, 1.69. The horizontal and
vertical differences of the major tune values between measurements and calculations are about 0.10 and 20.06, respectively. This discrepancy can be explained as follows.
The vertical tune has a large sensitivity to field gradient error because of the large beta function, as shown in Fig. 5.

<img src="./img/002-18.jpg"></img>

The focusing power of the fringe field, which is difficult
to calculate accurately, affects only the vertical tune.


A large excitation tune shift is found in the low excitation case shown in Fig. 15, and the tune shift at a high
excitation was relatively small. In the tracking in TOSCA’s
magnetic field, the tune shift is driven by the field deformation of iron saturation in the high magnetic field region.


Thus, this measured tune shift does not agree with the suggestion of the tune value from the TOSCA evaluation shown
in Fig. 7. Considering that the tune shift becomes large at
lower excitation currents, it seemed that this tune shift was
caused by an offset of the magnetic field. From simulations
based on the measured field under this offset assumption,
the expected offset of the magnetic field became about 5 G.
But this value was too large to be explained by either zero
offset of tesla meters or by the Earth’s magnetism. Thus,
it might be caused by the remanent magnetic field.


We tried to correct this tune shift with a one turn correction coil. This one turn correction coil was wound around
the magnet poles and its winding directions were opposite
between the F and D sectors, as shown in Fig. 16. Thus the
correction coil could change the relative excitation level of
the F sector against the D sector while keeping total bending power, because the total cross section of the correction
current loop of the F sector was almost equal to the cross
section of the correction current loop of the D sector. In
a first order approximation, the closed orbit deformation
was negligible and the tune shift was obtained from the
tune shift formula of Hill’s equation as


<img src="./img/002-19.jpg"></img>

where B0 and DB are the nominal value of the magnetic field in the designed orbit and the shift from the nominal value
at the F sector, respectively. The integrals R
F ds and R
D ds are the integration on the whole of the F and D sectors,
respectively. If the pole saturation is small, the ratio DBB0 is proportional to the ratio of the ampere turn product
between the correction coil and the main coil currents.

<img src="./img/002-20.jpg"></img>

The open box in Fig. 15 shows the operating point of
the magnetic field corrected by the correction coil with
a constant excitation current of 15 A. This correction
is effective for the low excitation level; then the evaluated
operating points come around the point 1.71, 1.69.
The used correction current of 15 A is within the designed
rating of the correction coil current of 100 A. Thus, 
by using the correction coil, we could suppress the 
excursion of the operating point sufficiently to use for the synchrotron.

### V. SUMMARY

We established an evaluation method for tune values
from the magnetic field distribution and developed a field
measurement scheme accurate enough to evaluate the tune
values. The evaluated operating point from the measured
field of the fabricated model magnet had reproducibility
within 0.005. The horizontal tune, which did not have a
large error propagation from the magnetic field as did the
vertical tune, agrees with the design value. In addition,
tune fluctuation with current excitation was found, but it
was confirmed that the tune fluctuations were easily 
compensated by a correction coil. Consequently, it was confirmed
that we could fabricate a combined function magnet
according to the design based on the three-dimensional
magnetic field calculations, by paying attention to parameters which have 
large error propagations. We cannot yet
confirm that a synchrotron ring works with this model magnet
because a complete magnet set to make such a synchrotron does not yet exist. However, from our tune value
evaluation of the measured magnetic field, we think that
the model magnet has sufficient performance for the construction of a synchrotron

### ACKNOWLEDGMENTS

略
