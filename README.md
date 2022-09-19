# Ritangle2021_Stage3
Solution for Stage 3 problem of Ritangle 2021 Puzzle competition

Competition home page: https://integralmaths.org/ritangle/2021
Puzzle home page: https://meiassets.blob.core.windows.net/integral-frontend-uploads/uploads/files/Ritangle_21_Stage_3_question.pdf

Puzzle statement:
-----------------
The final stage of Ritangle 2021 requires you to solve a similar problem but in three dimensions.
You have a cuboidal grid that measures u units × v units × w units where u, v and w are positive integers. Assume that the box occupies the space from (0, 0, 0) to (u, v, w).

You have to place a number of spheres in the grid, subject to the following rules:
1. Each sphere must be centred at one of the grid points. No two spheres may be centred on the same grid point.
2. The radius of each sphere must be a whole number of units.
3. Spheres may touch the boundary of the grid but may not cross it.
4. Spheres may touch one another but their boundaries may not cross. It is allowed for one sphere to be entirely inside another.
5. You may choose the values of 𝑢, 𝑣 and 𝑤, subject to the constraint 𝑢 + 𝑣 + 𝑤 ≤ 51.
6. You must maximise the total volume enclosed by all your spheres.

Express the total volume V as a multiple of the volume of a sphere with radius 1 unit. (i.e. calculate the volume in cubic units and divide by 4/3𝜋.).

My solution:
------------
I came up with 187 spheres located at various co-ordinates (see the file spheres.m or volumes.m for exact co-ordinates)

This matlab code checks my solution to the above problem. The puzzle is the last stage in a series of Maths puzzles competition called as 'Ritangle', hosted by Integral every year for Sixth Form students in England.

I used geogebra and other logic to nail down the exact co-ordinates and radii of the 187 spheres. I then used matlab code to check the sanctity of my solution, i.e. none of the spheres intersect any other sphere and they don't cross the boundaries.
