In our solution, we apply the formulas as they appear in the lectures.

We can see that the Euler method can makes our particles unstable (mainly when T is big) since, as T is bigger, the approximation of the differential
equation will be bigger. By setting a small T, the simulation is more stable, but more computation is necessary to have a fluid simulation.

We can also see that force-based walls are not efficient at all. This does not look like a real collision with a wall but more like some wind that would hinder particles
to escape the square.
Increasing the value of collision_stiffness solve the problem of particles going too far through the wall but makes our scene go unstable very easily.

Although we implemented the springs forces with Hooke's law and that the triangle shape is preserved, the shape in scene 5 is not necessarily preserved
when a strong force is applied to one of its node. Since no shape area preservation is implemented.

We needed around 4 hours to finish this lab (accounting the time to understand how to apply the lecture formulas)
The difficult part of this lab was understanding how the Euler method work and, as usual with graphical application, debugging when results are not correct.
