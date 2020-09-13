#version 100

#define MAX_MARCH_ITERATIONS 100
#define MAX_SCENE_DISTANCE 100.0
#define MIN_SURFACE_DISTANCE 0.01

precision mediump float;

uniform float u_time;
uniform vec2 u_resolution;

struct capsule_t {
    vec3 one_end;
    vec3 other_end;
    float radius;
};

float get_capsule_distance(vec3 point, capsule_t capsule) {
    vec3 point_to_top = point - capsule.one_end;
    vec3 capsule_axis = capsule.other_end - capsule.one_end;

    float t_along_axis = dot(point_to_top, capsule_axis) / dot(capsule_axis, capsule_axis);
    t_along_axis = clamp(t_along_axis, 0., 1.);
    vec3 point_projected_on_axis = capsule.one_end + capsule_axis*t_along_axis;

    return length(point - point_projected_on_axis) - capsule.radius;
}

// Returns the distance to the nearest object in a scene from a point.
//
// NOTA BENE: For now, we are hardcoding the two objects in the scene,
// a plane in the XZ plane and a sphere.
float get_distance(vec3 point) {
    // Note we store the radius of the sphere in the w coordinate.
    vec4 sphere = vec4(0., 1., 6., 1.);

    capsule_t capsule = capsule_t(vec3(-2., 1., 6.), vec3(-1., 2., 6.), 0.5);

    // The distance to the surface of a sphere is the
    // distance to its center minus its radius.
    float sphere_distance = length(point - sphere.xyz) - sphere.w;

    float capsule_distance = get_capsule_distance(point, capsule);

    // The distance is trivial to compute here.
    float plane_distance = point.y;

    // Choose the closest distance.
    return min(capsule_distance, plane_distance);
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
        float scene_distance = get_distance(new_point);

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
    float object_distance = get_distance(point);

    // This computes a normal by finding the distances between each of
    // three points slightly along the x, y, and z axes away from
    // the point. Each of these distances contributes to their
    // respective components of the normal vector, i.e.:
    //
    //                  dx*î + dy*ĵ + dz*k̂
    vec3 normal = object_distance - vec3(
        get_distance(point - vec3(0.01, 0., 0.)),
        get_distance(point - vec3(0., 0.01, 0.)),
        get_distance(point - vec3(0., 0., 0.01)));

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
 
    vec3 camera_position = vec3(0., 1., 0.);
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