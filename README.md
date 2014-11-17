swisscheese
===========

Computes the union, intersection, and subtraction areas for sets of circles.

I created this project to rapidly compute the area of masked apertures in astronomy for many different arbitrary arrangments of the masked areas. (A masked aperture, as I define it, is a circular area over which to sum the light accumulated on a detector, but with smaller circular areas masking out contaminating features removed.) However, this algorithm could have broader applications than astronomy, so I tried to make it as general as possible. As such, one can use it to take two arbitrary sets of circles, the union of each defnining its own shape, and compute the area of the intersection or difference of these sets, or the area of one with the other subtracted. (The area of the union may be computed by treating all of the circles as a single set, and this capability is available.)

----This is not optimized for speed with large numbers of circles.----
For my application, it was not necessary to create code that could handle large numbers of circles quickly. Rather, I needed accuracy and speed when dealing with lots of different configurations of just a few circles. The code meets those goals. For large numbers of circles, consider learning enough abstract geometry to implement the O(n log(n)) algorithm presented in 
http://www.sciencedirect.com/science/article/pii/0196677488900351#
