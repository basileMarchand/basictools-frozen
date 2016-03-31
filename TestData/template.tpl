****simulate
	***solver newton
		*iteration 30
		*ratio automatic 1.e-3
		%*divergence 2.
		%*security 1.2
		*algorithm p1p2p3
		%*automatic_time
		%*max_dtime 3.
		%*min_dtime 3.
		%*min_dtime 1.e-8

	***test {a}
		**load cycle 5
			*segment 200

time sig11 sig22 sig33 sig12 sig23 sig31
0.0  0.0   0.0   0.0   0.0   0.0   0.0
100. 250.  0.0   0.0   0.0   0.0   0.0
200. 0.0   0.0   0.0   0.0   0.0   0.0
300. -250.  0.0   0.0   0.0   0.0   0.0
400. 0.0  0.0   0.0   0.0   0.0   0.0

		**model
			*file {b}
			*integration runge_kutta 1.e-3 max 10000

		**no_binary_output
		**output
			*precision 8
			*dtime 1.
time sig11 eto11 evcum			
%time sig11 sig22 sig33 sig12 sig23 sig31 eto11 eto22 eto33 eto12 eto23 eto31 evcum X11 X22 X33 X12 X23 X31

****return