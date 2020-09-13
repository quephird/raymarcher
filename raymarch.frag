#version 100

#define MAX_MARCH_ITERATIONS 100
#define MAX_SCENE_DISTANCE 100.0
#define MIN_SURFACE_DISTANCE 0.01

precision mediump float;

uniform float u_time;
uniform vec2 u_resolution;

// Returns the distance to the nearest object in a scene from a point.
//
// NOTA BENE: For now, we are hardcoding the two objects in the scene,
// a plane in the XZ plane and a sphere.
float get_scene_distance(vec3 point) {
    // Note we store the radius of the sphere in the w coordinate.
    vec4 sphere = vec4(0., 1., 6., 1.);

    // The distance to the surface of a sphere is the
    // distance to its center minus its radius.
    float sphere_distance = length(point - sphere.xyz) - sphere.w;

    // The distance is trivial to compute here.
    float plane_distance = point.y;

    // Choose the closest distance.
    return min(sphere_distance, plane_distance);
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
        float scene_distance = get_scene_distance(new_point);

        // Add that distance to the current one
        ray_distance += scene_distance;

        // If we've gone too far or we're sufficiently close to an object
        // stop iterating.
        if(ray_distance > MAX_SCENE_DISTANCE || scene_distance < MIN_SURFACE_DISTANCE)
            break;
    }
    
    return ray_distance;
}

void main() {
    // Center image and set aspect ratio to something pleasing
    vec2 uv = (gl_FragCoord.xy - 0.5*u_resolution) / u_resolution.y;

    vec3 color = vec3(0.0);
 
    vec3 camera_position = vec3(0., 1., 0.);
    vec3 camera_direction = normalize(vec3(uv.xy, 1.0));

    // See if we get a hit to an object
    float object_distance = march(camera_position, camera_direction);

    // Since all objects in this particular scene at at least one unit away,
    // we need to attenuate the magnitude of `object_distance` if we
    // want to use it as the color.
    color = vec3(object_distance/6.);

    // Set pixel color!!!
    gl_FragColor = vec4(color, 1.0);
}