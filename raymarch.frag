#version 100

#define MAX_MARCH_ITERATIONS 100
#define MAX_SCENE_DISTANCE 100.0
#define MIN_SURFACE_DISTANCE 0.01

#define RED vec3(1.0, 0.0, 0.0)
#define ORANGE vec3(1.0, 0.5, 0.0)
#define YELLOW vec3(1.0, 1.0, 0.0)
#define GREEN vec3(0.0, 1.0, 0.0)
#define BLUE vec3(0.0, 0.0, 1.0)
#define PURPLE vec3(0.5, 0.0, 1.0)

precision mediump float;

uniform float u_time;
uniform vec2 u_resolution;

struct nearest_object_t {
    float d;
    vec3 color;
};

struct rounded_box_t {
    vec3 position;
    vec3 dimensions;
    float outer_radius;
    vec3 color;
};

//
//                                 _─‾‾|‾‾‾‾‾‾‾‾‾|‾‾─_
//                               _‾  r |         |    ‾_
//                              |______|_________|______|
//                              |      |         |      |
//                              |      |         |      |
//                              |      | w       |      |
//                              |      |────· B  |      |
//                              |      |    |    |      |
//                              |      |   h|    |      |
//                              |      |    |    |      |
//                              |‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾|
//                               ‾_    |         |    _‾
//                                 ‾─__|_________|__─‾
//

nearest_object_t get_rounded_box_distance(vec3 point, rounded_box_t rounded_box) {
    vec3 half_sizes = rounded_box.dimensions / 2.0;
    vec3 q = abs(point - rounded_box.position) - half_sizes;
    float d = length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - rounded_box.outer_radius;

    return nearest_object_t(d, rounded_box.color);
}

struct box_t {
    vec3 position;
    vec3 dimensions;
    vec3 color;
};


//
//   2)                         |‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾|
//                              |                       |
//                              |                       |
//                              |    w/2                |
//                              |───────────· B         |
//                              |           |           |
//                              |           | h/2       |
//                              |           |           |
//                              *___________|___________|
//                     d ___───‾|
//              ___───‾‾‾       | dy
//        *──‾─‾────────────────*
//     P              dx
//
//   3)                         |‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾|
//                              |                       |
//                              |                       |
//                              |    w/2                |
//                              |───────────· B         |
//     P           d = dx       |           |           |
//     *─────_─_────────────────*           | h/2       |
//              ‾‾‾───___       | dy        |           |
//                       ‾‾‾────|___________|___________|
//
//
//   4)                         |‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾|
//                              |                       |
//                              |                       |
//                              |    w/2                |
//                              |───────────· B         |
//                              |           |           |
//                              |           | h/2       |
//                              |  dx       |           |
//                              |______*____|___________|
//                               \     |
//                                \    |
//                                 \   | d = dy
//                                  \  |
//                                   \ |
//                                    \|
//                                     * P
//
//  Let's just consider a square for a moment. The minimum distance d from P
//  to the square centered at B will depend on the two distances, dx and dy.
//  There are four possibilities, three of which are pictured above:
//
//  1) dx<0 and dy<0: point is inside the cube, which we will ignore for now
//  2) dx<0 and dy>0: distance to square is dy
//  3) dx>0 and dy<0: distance to square is dx
//  4) dx>0 and dy>0: distance to square is d = √(dx² + dy²)
//
//  We can express this very pithily by noting that in each case, taking the
//  `max` of the [dx, dy] and [0, 0] vectors will result in a vector whose
//  components are nonnegative. (This illustrates how `max` can be exploited
//  to avoid `if/else` conditionals.) That is:
//
//  2) dx<0 and dy>0: max([dx, dy], [0, 0]) -> [dx, 0]
//  3) dx>0 and dy<0: max([dx, dy], [0, 0]) -> [0, dy]
//  4) dx>0 and dy>0: max([dx, dy], [0, 0]) -> [dx, dy]
//
//  In each case, we can take the length of the resulting vector to get the
//  desired distance, either to a corner or to a face. It's the same thing
//  when we move to three dimensions.

nearest_object_t get_box_distance(vec3 point, box_t box) {
    vec3 half_sizes = box.dimensions / 2.0;
    float d = length(max(abs(point - box.position) - half_sizes, vec3(0.0)));

    return nearest_object_t(d, box.color);
}

struct torus_t {
    vec3 position;
    float primary_radius;
    float secondary_radius;
    vec3 color;
};

//
//                        __─────‾‾‾‾‾─────__
//                    _─‾‾                   ‾‾─_
//       P           /  _─_                      \
//       *_         |  | | |      _* T            |
//       | ‾‾──__    \│  |r │ _──‾               /
//       |       ‾‾──_:_ | _│‾   R            _─‾
//       |           D│_─‾‾─:────_____──────‾‾
//       |         _──|  C  |
//       |     _──‾    |   |
//       | _──‾         ‾─‾
//       *‾
//       P'
//
//  So, this time, the minimal distance to a torus lying in the xz plane, the length
//  of line segment PD, requires a few intermediate calculations. The first thing we
//  can calculate is the length of the segment from the point on the xz plane that P
//  projects onto, namely P', to the center of the torus, T. That means that the length
//  of P' to the center of the closest outer circle, namely C, is the length of P'T
//  minus the outer radius R. Since PP' abd P'C form a right triangle, and the length
//  of PP' is simply the y value of P, then the length of PC can be easily calculated
//  via the Pythagorean theorem. subtracting the outer radius, r, yields:
//
//                 length(PD) = √(length(P'C)^2 + P.y^2)) - r
//

nearest_object_t get_torus_distance(vec3 point, torus_t torus) {
    float dx = length(point.xz - torus.position.xz) - torus.primary_radius;
    float dy = point.y - torus.position.y;

    float d = sqrt(dx*dx + dy*dy) - torus.secondary_radius;
    return nearest_object_t(d, torus.color);
}

struct cylinder_t {
    vec3 one_end;
    vec3 other_end;
    float radius;
    vec3 color;
};

//     1)                      _───‾‾‾───_
//                            |_    * A  _|
//                            | ‾‾‾───‾‾‾ |
//                            |     |     |
//                            |     |     |
//                            |     |     |
//                            |     |     |
//                            |  r  |     |
//                            |_____* B  _|
//               d      __──‾‾| ‾‾‾───‾‾‾
//                __──‾‾      |     |
//          __──‾‾            |     |
//     *──‾‾──────────────────|─────|
//     P
//
//     2)                      _───‾‾‾───_
//                            |_    * A  _|
//                            | ‾‾‾───‾‾‾ |
//                 d          |  r  |     |
//     *──────────────────────|─────|     |
//     P                      |     |     |
//                            |     |     |
//                            |     |     |
//                            |_    * B  _|
//                              ‾‾‾───‾‾‾
//
//
//     3)                      _───‾‾‾───_
//                            |_    * A  _|
//                            | ‾‾‾───‾‾‾ |
//                            |     |     |
//                            |     |     |
//                            |     |     |
//                            |     |     |
//                            |  r  |     |
//                            |_____* B  _|
//                             \‾‾‾───‾‾‾
//                              |  |
//                              |  |
//                               \ | d
//                                ||
//                                 |
//                               P *
//
//  There are four possibilities, three of which are pictured above:
//
//  1) point and projection of point are external _infinite_ cylinder
//  2) point is external to but projection of point is inside cylinder
//  3) point is inside _infinite_ cylinder
//  4) point is inside the cylinder, which we will ignore for now
//

nearest_object_t get_cylinder_distance(vec3 point, cylinder_t cylinder) {
    vec3 point_to_top = point - cylinder.one_end;
    vec3 cylinder_axis = cylinder.other_end - cylinder.one_end;
    float axis_length = length(cylinder_axis);

    float t_along_axis = dot(point_to_top, cylinder_axis) / dot(cylinder_axis, cylinder_axis);
    vec3 point_projected_on_axis = cylinder.one_end + cylinder_axis*t_along_axis;

    float distance_to_infinite_cylinder = length(point - point_projected_on_axis) - cylinder.radius;
    float distance_to_cylinder_along_axis = (abs(t_along_axis - 0.5) - 0.5) * axis_length;
    float distance_to_cylinder_edge = length(max(vec2(distance_to_infinite_cylinder, distance_to_cylinder_along_axis), 0.0));
    float internal_distance = min(max(distance_to_infinite_cylinder, distance_to_cylinder_along_axis), 0.0);

    float d = internal_distance + distance_to_cylinder_edge;

    return nearest_object_t(d, cylinder.color);
}

struct capsule_t {
    vec3 one_end;
    vec3 other_end;
    float radius;
    vec3 color;
};

//
//                                   ,──‾‾‾──,
//                                 ,'         `.
//                               __|_────* A   |
//                _____─────‾‾‾‾‾  |     |     |
//     *─────‾─‾─‾─────────────────|─────* C   |
//     P                           |     |     |
//                                 |─────* B   |
//                                 `. r       ,'
//                                   `──___──'
//
//  When P is such that its perpendicular falls between both capsule centers,
//  A and B, intersecting at C, then the value of t along AB is computed
//  by projecting PA onto AB, in which case t will be in the interval [0, 1].
//
//     *───__────────────────────────────* C
//     P     ‾‾──__                      ·
//                 ‾‾──__                ·
//                       ‾‾──__      ,──‾·‾──,
//                             ‾‾──,'_   ·    `.
//                                 |  ‾‾─* A   |
//                                 |     |     |
//                                 |     |     |
//                                 |     |     |
//                                 |─────* B   |
//                                 `. r       ,'
//                                   `──___──'
//
// However when PC does _not_ intersect between A and B, then t falls _outside_
// [0, 1], in the example above t < 0, and so we need to ensure we clamp t to [0, 1].
// Doing that will effectively choose the closer of the two centers, and the
// distance formula still holds.
//
nearest_object_t get_capsule_distance(vec3 point, capsule_t capsule) {
    vec3 point_to_top = point - capsule.one_end;
    vec3 capsule_axis = capsule.other_end - capsule.one_end;

    float t_along_axis = dot(point_to_top, capsule_axis) / dot(capsule_axis, capsule_axis);
    t_along_axis = clamp(t_along_axis, 0., 1.);
    vec3 point_projected_on_axis = capsule.one_end + capsule_axis*t_along_axis;

    float d = length(point - point_projected_on_axis) - capsule.radius;

    return nearest_object_t(d, capsule.color);
}

struct sphere_t {
    vec3 position;
    float radius;
    vec3 color;
};

//
//     *─────_____                   _──‾|‾──_
//     P          ‾‾‾‾‾─────_____   /    |r   \
//                               ‾‾|‾────*     |
//                                  \    S    /
//                                   ‾──___──‾
//
//  The distance between the point P and the sphere centered at S
//  is the length of PS minus the radius r.
//
nearest_object_t get_sphere_distance(vec3 point, sphere_t sphere) {
    float d = length(point - sphere.position) - sphere.radius;

    return nearest_object_t(d, sphere.color);
}

struct scene_t {
    sphere_t sphere;
    capsule_t capsule;
    torus_t torus;
    box_t box;
    cylinder_t cylinder;
    rounded_box_t rounded_box;
};

nearest_object_t closer_object(nearest_object_t object1, nearest_object_t object2) {
    if (object1.d < object2.d) {
        return object1;
    } else {
        return object2;
    }
}

// Returns the distance to the nearest object in a scene from a point.
//
// NOTA BENE: For now, we are hardcoding some knowledge of the scene here.
nearest_object_t get_nearest_object(vec3 point, scene_t scene) {
    nearest_object_t sphere = get_sphere_distance(point, scene.sphere);
    nearest_object_t capsule = get_capsule_distance(point, scene.capsule);
    nearest_object_t torus= get_torus_distance(point, scene.torus);
    nearest_object_t box = get_box_distance(point, scene.box);
    nearest_object_t cylinder = get_cylinder_distance(point, scene.cylinder);
    nearest_object_t rounded_box = get_rounded_box_distance(point, scene.rounded_box);

    // TODO: promote plane to first class object
    // The distance is trivial to compute here.
    float plane_distance = point.y;
    nearest_object_t plane = nearest_object_t(plane_distance, vec3(1.0, 1.0, 1.0));

    return closer_object(plane,
               closer_object(capsule,
                   closer_object(box,
                       closer_object(torus,
                           closer_object(cylinder,
                               closer_object(rounded_box, sphere))))));
}

// This takes a ray, whose origin and direction are passed in,
// and returns the distance to the closest object in a scene.
nearest_object_t march(vec3 ray_origin, vec3 ray_direction, scene_t scene) {
    // Start at the ray origin...
	float ray_distance = 0.;
    vec3 object_color = vec3(0.0);

    for(int i=0; i<MAX_MARCH_ITERATIONS; i++) {
        // March down the ray the current distance
    	vec3 new_point = ray_origin + ray_direction*ray_distance;

        // Get the distance to the closest object
        nearest_object_t scene_distance = get_nearest_object(new_point, scene);

        // Add that distance to the current one
        ray_distance += scene_distance.d;
        object_color = scene_distance.color;

        // If we've gone too far or we're sufficiently close to an object
        // stop iterating.
        if(ray_distance > MAX_SCENE_DISTANCE || scene_distance.d < MIN_SURFACE_DISTANCE)
            break;
    }

    return nearest_object_t(ray_distance, object_color);
}

// This function computes the normal vector at the point passed in.
// NOTA BENE: `get_distance` is what's doing most of the work here
// since only it knows about the objects in the scene.
vec3 get_normal_vector(vec3 point, scene_t scene) {
    nearest_object_t object_distance = get_nearest_object(point, scene);

    // This computes a normal by finding the distances between each of
    // three points slightly along the x, y, and z axes away from
    // the point. Each of these distances contributes to their
    // respective components of the normal vector, i.e.:
    //
    //                  dx*î + dy*ĵ + dz*k̂
    nearest_object_t nearest_x = get_nearest_object(point - vec3(0.01, 0., 0.), scene);
    nearest_object_t nearest_y = get_nearest_object(point - vec3(0., 0.01, 0.), scene);
    nearest_object_t nearest_z = get_nearest_object(point - vec3(0., 0., 0.01), scene);

    vec3 normal = object_distance.d - vec3(nearest_x.d, nearest_y.d, nearest_z.d);

    return normalize(normal);
}

// This function returns a scalar float representing
// the brightness of color at the point passed in, given
// a single white light source that is hardcoded below.
//
// Note that the light is brightest directly above a surface point,
// and darker as the angle between the two vectors increases.
//
//                              * light source
//
//         normal vector
//                       |   /
//                       |  /  light vector
//                       | /
//                       |/
//           ____________*______________
//                  point
//
float get_diffused_light(vec3 point, scene_t scene) {
    vec3 diffused_color = vec3(0.0);

    // Single light source
    vec3 light_position = vec3(4.0*cos(u_time), 5., 3.0+4.0*sin(u_time));
    vec3 light_direction = normalize(light_position - point);

    // For now, just compute a color based on the magnitude
    // of the dot product
    vec3 normal = get_normal_vector(point, scene);
    float color = dot(normal, light_direction);

    // Now we check to see if another object is between the given
    // point and the light source.
    nearest_object_t object_distance = march(point+normal*MIN_SURFACE_DISTANCE, light_direction, scene);

    // If so, then darken the color to simulate shadowing
    if (object_distance.d < length(light_position - point))
        color *= 0.1;

    // The dot product above actually returns a value in [-1, 1],
    // and so we need to make sure the returned value is in [0, 1].
    return clamp(color, 0., 1.);
}

void main() {
    // Center image and set aspect ratio to something pleasing
    vec2 uv = (gl_FragCoord.xy - 0.5*u_resolution) / u_resolution.y;

    vec3 color = vec3(0.0);

    // Camera position rotates about y-axis now
    vec3 camera_position = vec3(-10.*cos(u_time+1.57), 3.0, 4.-10.*sin(u_time+1.57));

    // Camera direction now always faces the center of the scene!
    mat3 rotation_about_y_axis;
    rotation_about_y_axis[0] = vec3(cos(u_time), 0.0, sin(u_time));
    rotation_about_y_axis[1] = vec3(0.0, 1.0, 0.0);
    rotation_about_y_axis[2] = vec3(-sin(u_time), 0.0, cos(u_time));
    vec3 camera_direction = normalize(rotation_about_y_axis*vec3(uv.xy, 1.0));

    // Construct the scene here
    sphere_t sphere = sphere_t(vec3(1.5, 1.5, 6.), 1.0, RED);
    capsule_t capsule = capsule_t(vec3(-2., 1., 6.), vec3(-1., 2., 6.), 0.5, GREEN);
    torus_t torus = torus_t(vec3(0., 0.5, 4.), 1., 0.25, ORANGE);
    box_t box = box_t(vec3(-2.5, 1., 4.), vec3(1.0, 1.0, 1.0), PURPLE);
    cylinder_t cylinder = cylinder_t(vec3(3., 0.5, 4.), vec3(2., 1.5, 3.), 0.5, YELLOW);
    rounded_box_t rounded_box = rounded_box_t(vec3(1.5, 0.7, 2.), vec3(1., 1., 1.), 0.1, BLUE);

    scene_t scene = scene_t(sphere, capsule, torus, box, cylinder, rounded_box);

    // See if we get a hit to an object
    nearest_object_t object_distance = march(camera_position, camera_direction, scene);

    // Travel down that distance to the point
    vec3 surface_point = camera_position + camera_direction*object_distance.d;

    // Compute grey-scale color
    float diffused_light = get_diffused_light(surface_point, scene);
    color = object_distance.color*0.5 + diffused_light;

    // Set pixel color!!!
    gl_FragColor = vec4(color, 1.0);
}