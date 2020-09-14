#version 100

#define MAX_MARCH_ITERATIONS 100
#define MAX_SCENE_DISTANCE 100.0
#define MIN_SURFACE_DISTANCE 0.01

precision mediump float;

uniform float u_time;
uniform vec2 u_resolution;

struct torus_t {
    vec3 position;
    float primary_radius;
    float secondary_radius;
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
float get_torus_distance(vec3 point, torus_t torus) {
    float dx = length(point.xz - torus.position.xz) - torus.primary_radius;
    float dy = point.y - torus.position.y;

    return sqrt(dx*dx + dy*dy) - torus.secondary_radius;
}

struct capsule_t {
    vec3 one_end;
    vec3 other_end;
    float radius;
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
float get_capsule_distance(vec3 point, capsule_t capsule) {
    vec3 point_to_top = point - capsule.one_end;
    vec3 capsule_axis = capsule.other_end - capsule.one_end;

    float t_along_axis = dot(point_to_top, capsule_axis) / dot(capsule_axis, capsule_axis);
    t_along_axis = clamp(t_along_axis, 0., 1.);
    vec3 point_projected_on_axis = capsule.one_end + capsule_axis*t_along_axis;

    return length(point - point_projected_on_axis) - capsule.radius;
}

struct sphere_t {
    vec3 position;
    float radius;
};

//
//     *─────_____                   ,──‾|‾──,
//     P          ‾‾‾‾‾─────_____  ,'    |r   `.
//                               ‾‾|‾────*     |
//                                 `.    S    ,'
//                                   `──___──'
//
//  The distance between the point P and the sphere centered at S
//  is the length of PS minus the radius r.
//
float get_sphere_distance(vec3 point, sphere_t sphere) {
    return length(point - sphere.position) - sphere.radius;
}

// Returns the distance to the nearest object in a scene from a point.
//
// NOTA BENE: For now, we are hardcoding objects in the scene
float get_nearest_distance(vec3 point) {
    sphere_t sphere = sphere_t(vec3(1.5, 1.5, 6.), 1.0);
    capsule_t capsule = capsule_t(vec3(-2., 1., 6.), vec3(-1., 2., 6.), 0.5);
    torus_t torus = torus_t(vec3(0., 0.5, 4.), 1., 0.25);

    float sphere_distance = get_sphere_distance(point, sphere);
    float capsule_distance = get_capsule_distance(point, capsule);
    float torus_distance = get_torus_distance(point, torus);

    // The distance is trivial to compute here.
    float plane_distance = point.y;

    // Choose the closest distance.
    return min(plane_distance,
               min(sphere_distance,
                   min(capsule_distance, torus_distance)));
}

// This takes a ray, whose origin and direction are passed in,
// and returns the distance to the closest object in a scene.
float march(vec3 ray_origin, vec3 ray_direction) {
    // Start at the ray origin...
	float ray_distance = 0.;

    for(int i=0; i<MAX_MARCH_ITERATIONS; i++) {
        // March down the ray the current distance
    	vec3 new_point = ray_origin + ray_direction*ray_distance;

        // Get the distance to the closest object
        float scene_distance = get_nearest_distance(new_point);

        // Add that distance to the current one
        ray_distance += scene_distance;

        // If we've gone too far or we're sufficiently close to an object
        // stop iterating.
        if(ray_distance > MAX_SCENE_DISTANCE || scene_distance < MIN_SURFACE_DISTANCE)
            break;
    }

    return ray_distance;
}

// This function computes the normal vector at the point passed in.
// NOTA BENE: `get_distance` is what's doing most of the work here
// since only it knows about the objects in the scene.
vec3 get_normal_vector(vec3 point) {
    float object_distance = get_nearest_distance(point);

    // This computes a normal by finding the distances between each of
    // three points slightly along the x, y, and z axes away from
    // the point. Each of these distances contributes to their
    // respective components of the normal vector, i.e.:
    //
    //                  dx*î + dy*ĵ + dz*k̂
    vec3 normal = object_distance - vec3(
        get_nearest_distance(point - vec3(0.01, 0., 0.)),
        get_nearest_distance(point - vec3(0., 0.01, 0.)),
        get_nearest_distance(point - vec3(0., 0., 0.01)));

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
float get_diffused_light(vec3 point) {
    // Single light source
    vec3 light_position = vec3(0., 5., 6.);
    vec3 light_direction = normalize(light_position - point);

    // For now, just compute a color based on the magnitude
    // of the dot product
    vec3 normal = get_normal_vector(point);
    float color = dot(normal, light_direction);

    // Now we check to see if another object is between the given
    // point and the light source.
    float object_distance = march(point+normal*MIN_SURFACE_DISTANCE, light_direction);

    // If so, then darken the color to simulate shadowing
    if (object_distance < length(light_position - point))
        color *= 0.1;

    // The dot product above actually returns a value in [-1, 1],
    // and so we need to make sure the returned value is in [0, 1].
    return clamp(color, 0., 1.);
}

void main() {
    // Center image and set aspect ratio to something pleasing
    vec2 uv = (gl_FragCoord.xy - 0.5*u_resolution) / u_resolution.y;

    vec3 color = vec3(0.0);

    vec3 camera_position = vec3(0., 2., -2.);
    vec3 camera_direction = normalize(vec3(uv.xy, 1.0));

    // See if we get a hit to an object
    float object_distance = march(camera_position, camera_direction);

    // Travel down that distance to the point
    vec3 surface_point = camera_position + camera_direction*object_distance;

    // Compute grey-scale color
    float diffused_light = get_diffused_light(surface_point);
    color = vec3(diffused_light);

    // Set pixel color!!!
    gl_FragColor = vec4(color, 1.0);
}