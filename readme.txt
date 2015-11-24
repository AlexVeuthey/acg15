In our solution, we apply the force calculation for the lecture week9 and the
impulse based collisons from the slides of week10.

We saw that in the rigid body class, the points vector of vec2 correspond to the 
points of the shape after transformation, as the r, vector of relative positions is
the relative positions BEFORE the transformation. This has caused some confusions
about our intermediate results. The fact that r is the original value can be seen in
update_points() method of rigid_body where "points" is modified but not "r".

In our impulse based collisions we also have set the "elasticity" value to 0.9. This has been
done to have the best visual effect, an elasticity of 1.0 make the shape bounce too much
and an elasticity of near 0.0 is dull as the shape doesn't bounce.

In our implementation of the impulse based collision, we took care of checking the
value of v_rel, the relative speed of the particle in relation to the wall. If this
value is smaller than 0, then the particle is moving inside the wall in the wrong direction
and thus, we apply the change of velocity. Otherwise the particle is going in the right
direction (or is still on the line) and we do nothing. This allows us to have shapes that
stay still on the ground.
