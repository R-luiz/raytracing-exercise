
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    double x, y, z;
} vec3;

// Vector addition
vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){ a.x + b.x, a.y + b.y, a.z + b.z };
}

// Vector subtraction
vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){ a.x - b.x, a.y - b.y, a.z - b.z };
}

// Dot product
double vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Multiply vector by scalar
vec3 vec3_scale(vec3 v, double s) {
    return (vec3){ v.x * s, v.y * s, v.z * s };
}

// Vector length squared
double vec3_length_squared(vec3 v) {
    return vec3_dot(v, v);
}

typedef struct {
    vec3 origin;
    vec3 direction;
} ray;

// Compute a point along a ray
vec3 ray_at(ray r, double t) {
    return vec3_add(r.origin, vec3_scale(r.direction, t));
}

// Tests if a ray hits a sphere.
/*
double hit_sphere(vec3 sphere_center, double radius, ray r) {
    vec3 oc = vec3_sub(r.origin, sphere_center);
    double a = vec3_length_squared(r.direction);
    double half_b = vec3_dot(oc, r.direction);
    double c = vec3_length_squared(oc) - radius * radius;
    double discriminant = half_b * half_b - a * c;

    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant)) / a;
    }
}
*/
vec3 vec3_unit(vec3 v) {
    double length = sqrt(vec3_length_squared(v));
    return vec3_scale(v, 1.0 / length);
}

// Computes the color of a ray.
/*
vec3 ray_color(ray r) {
    double t = hit_sphere((vec3){0, 0, -1}, 0.5, r);
    if (t > 0.0) {
        vec3 N = vec3_unit(vec3_sub(ray_at(r, t), (vec3){0, 0, -1}));
        return vec3_scale((vec3){N.x + 1, N.y + 1, N.z + 1}, 0.5);
    }

    vec3 unit_direction = vec3_unit(r.direction);
    t = 0.5 * (unit_direction.y + 1.0);
    return vec3_add(vec3_scale((vec3){1.0, 1.0, 1.0}, 1.0 - t), vec3_scale((vec3){0.5, 0.7, 1.0}, t));
}*/

typedef struct {
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
} camera;

camera camera_init(double aspect_ratio) {
    camera cam;
    double viewport_height = 2.0;
    double viewport_width = aspect_ratio * viewport_height;
    double focal_length = 1.0;

    cam.origin = (vec3){0, 0, 0};
    cam.horizontal = (vec3){viewport_width, 0, 0};
    cam.vertical = (vec3){0, viewport_height, 0};
    cam.lower_left_corner = vec3_sub(
        vec3_sub(
            vec3_sub(cam.origin, vec3_scale(cam.horizontal, 0.5)),
            vec3_scale(cam.vertical, 0.5)
        ),
        (vec3){0, 0, focal_length}
    );

    return cam;
}

ray get_ray(camera cam, double u, double v) {
    return (ray){
        cam.origin,
        vec3_sub(
            vec3_add(
                vec3_add(cam.lower_left_corner, vec3_scale(cam.horizontal, u)),
                vec3_scale(cam.vertical, v)
            ),
            cam.origin
        )
    };
}

typedef struct {
    vec3 p;
    vec3 normal;
    double t;
} hit_record;

typedef struct {
    vec3 center;
    double radius;
} Sphere;

typedef struct {
    Sphere *spheres;
    int size;
} hittable_list;


int hit_sphere(const Sphere *sphere, const ray *r, double t_min, double t_max, hit_record *rec) {
    vec3 oc = vec3_sub(r->origin, sphere->center);
    double a = vec3_length_squared(r->direction);
    double half_b = vec3_dot(oc, r->direction);
    double c = vec3_length_squared(oc) - sphere->radius * sphere->radius;
    double discriminant = half_b * half_b - a * c;

    if (discriminant > 0) {
        double sqrtd = sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        double root = (-half_b - sqrtd) / a;
        if (root < t_min || root > t_max) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || root > t_max) {
                return 0;
            }
        }

        rec->t = root;
        rec->p = ray_at(*r, rec->t);
        rec->normal = vec3_scale(vec3_sub(rec->p, sphere->center), 1.0 / sphere->radius);
        return 1;
    }

    return 0;
}

int hit(const hittable_list *world, const ray *r, double t_min, double t_max, hit_record *rec) {
    hit_record temp_rec;
    int hit_anything = 0;
    double closest_so_far = t_max;

    for (int i = 0; i < world->size; i++) {
        if (hit_sphere(&world->spheres[i], r, t_min, closest_so_far, &temp_rec)) {
            hit_anything = 1;
            closest_so_far = temp_rec.t;
            *rec = temp_rec;
        }
    }

    return hit_anything;
}

vec3 ray_color(const ray *r, const hittable_list *world) {
    hit_record rec;
    if (hit(world, r, 0, INFINITY, &rec)) {
        return vec3_scale((vec3){rec.normal.x + 1, rec.normal.y + 1, rec.normal.z + 1}, 0.5);
    } else {
        // Background color
        vec3 unit_direction = vec3_unit(r->direction);
        double t = 0.5 * (unit_direction.y + 1.0);
        return vec3_add(vec3_scale((vec3){1.0, 1.0, 1.0}, 1.0 - t), vec3_scale((vec3){0.5, 0.7, 1.0}, t));
    }
}

int main() {
    // Image
    const double aspect_ratio = 16.0 / 9.0;
    const int image_width = 800;
    const int image_height = (int)(image_width / aspect_ratio);
    char name[20];

    // Camera
    camera cam = camera_init(aspect_ratio);

    // File handling
    srand(time(NULL));
    sprintf(name, "image_%d.ppm", rand() % 9999);
    FILE *file = fopen(name, "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Write the header to the file
    fprintf(file, "P3\n%d %d\n255\n", image_width, image_height);

    // World
    hittable_list world;
    Sphere spheres[3];
    spheres[0] = (Sphere){.center = {0, 0, -1}, .radius = 0.5};
    spheres[1] = (Sphere){.center = {0, -100.5, -1}, .radius = 100};
    spheres[2] = (Sphere){.center = {1, 0, -1}, .radius = 0.4};
    world.spheres = spheres;
    world.size = 3;

    // Render and write to file
    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            double u = (double)i / (image_width - 1);
            double v = (double)j / (image_height - 1);
            ray r = get_ray(cam, u, v);
            vec3 pixel_color = ray_color(&r, &world);
            fprintf(file, "%d %d %d\n",
                    (int)(255.999 * pixel_color.x),
                    (int)(255.999 * pixel_color.y),
                    (int)(255.999 * pixel_color.z));
        }
    }

    // Close the file
    fclose(file);
    return 0;
}

