/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: rluiz <rluiz@student.42.fr>                +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2023/12/23 21:24:56 by rluiz             #+#    #+#             */
/*   Updated: 2023/12/24 02:24:31 by rluiz            ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct s_vec3
{
	double		x;
	double		y;
	double		z;
}				t_vec3;

typedef t_vec3	t_point3;
typedef struct s_ray
{
	t_point3	origin;
	t_vec3		direction;
}				t_ray;

t_vec3	vec3_add(t_vec3 v, t_vec3 w)
{
	t_vec3	ret;

	ret = (t_vec3){.x = v.x + w.x, .y = v.y + w.y, .z = v.z + w.z};
	return (ret);
}

t_vec3	vec3_sub(t_vec3 v, t_vec3 w)
{
	t_vec3	ret;

	ret = (t_vec3){.x = v.x - w.x, .y = v.y - w.y, .z = v.z - w.z};
	return (ret);
}

t_vec3	vec3_scale(t_vec3 v, double scalar)
{
	t_vec3	ret;

	ret = (t_vec3){.x = v.x * scalar, .y = v.y * scalar, .z = v.z * scalar};
	return (ret);
}

double	vec3_magnitude(t_vec3 v)
{
	return (sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

t_vec3	vec3_unit(t_vec3 v)
{
	t_vec3	ret;
	double	m;

	m = vec3_magnitude(v);
	ret = (t_vec3){.x = v.x / m, .y = v.y / m, .z = v.z / m};
	return (ret);
}

double	vec3_dot(t_vec3 v, t_vec3 w)
{
	return (v.x * w.x + v.y * w.y + v.z * w.z);
}

t_vec3	vec3_cross(t_vec3 v, t_vec3 w)
{
	t_vec3	ret;

	ret = (t_vec3){.x = v.y * w.z - v.z * w.y, .y = v.z * w.x - v.x * w.z,
		.z = v.x * w.y - v.y * w.x};
	return (ret);
}

t_point3	ray_at(t_ray ray, double t)
{
	return (vec3_add(ray.origin, vec3_scale(ray.direction, t)));
}

void	vec3_print(t_vec3 v, FILE *file)
{
	fprintf(file, "%d %d %d\n", (int)(v.x * 255.99), (int)(v.y * 255.99),
		(int)(v.z * 255.99));
}

double	hit_sphere(t_point3 *center, double radius, t_ray *ray)
{
	t_vec3	oc;
	double	a;
	double	b;
	double	c;
	double	discriminant;

	oc = vec3_sub(ray->origin, *center);
	a = vec3_dot(ray->direction, ray->direction);
	b = 2.0 * vec3_dot(oc, ray->direction);
	c = vec3_dot(oc, oc) - radius * radius;
	discriminant = b * b - 4 * a * c;
	if (discriminant < 0)
		return (-1);
	else
		return (-b - sqrt(discriminant) / (2.0 * a));
}

t_vec3	ray_color(t_ray *ray)
{
	t_vec3	unit_direction;
	double	t;
	double	a;

	t = hit_sphere(&(t_point3){.x = 0.5, .y = 0.2, .z = -1}, 0.5, ray);
	if (t > 0.0)
	{
		unit_direction = vec3_unit(vec3_sub(ray_at(*ray, t), (t_point3){.x = 0,
				.y = 0, .z = -1}));
		return (vec3_scale((t_vec3){.x = unit_direction.x + 1.0,
				.y = unit_direction.y + 1.0, .z = unit_direction.z + 1.0},
				0.5));
	}
	unit_direction = vec3_unit(ray->direction);
	a = 0.5 * (unit_direction.y + 1.0);
	return (vec3_add(vec3_scale((t_vec3){.x = 1.0, .y = 1.0, .z = 1.0}, 1.0
				- a), vec3_scale((t_vec3){.x = 0.2, .y = 0.4, .z = 1.0}, a)));
}

int	main(void)
{
	int			image_width;
	int			image_height;
	char		*name;
	FILE		*file;
	t_vec3		red;
	t_vec3		green;
	t_vec3		blue;
	t_vec3		white;
	t_vec3		black;
	t_vec3		yellow;
	t_vec3		magenta;
	t_vec3		cyan;
	t_vec3		color;
	double		r;
	double		g;
	double		b;
	double		aspect_ratio;
	int			viewport_height;
	int			viewport_width;
	double		focal_length;
	t_vec3		camera_center;
	t_vec3		viewport_u;
	t_vec3		viewport_v;
	t_vec3		pixel_delta_u;
	t_vec3		pixel_delta_v;
	t_vec3		viewport_upper_left;
	t_vec3		pixel00;
	t_point3	pixel_center;
	t_vec3		ray_direction;
	t_ray		ray;

	red = (t_vec3){.x = 1, .y = 0, .z = 0};
	green = (t_vec3){.x = 0, .y = 1, .z = 0};
	blue = (t_vec3){.x = 0, .y = 0, .z = 1};
	white = (t_vec3){.x = 1, .y = 1, .z = 1};
	black = (t_vec3){.x = 0, .y = 0, .z = 0};
	yellow = (t_vec3){.x = 1, .y = 1, .z = 0};
	magenta = (t_vec3){.x = 1, .y = 0, .z = 1};
	cyan = (t_vec3){.x = 0, .y = 1, .z = 1};
	// Image
	image_width = 400;
	aspect_ratio = 16.0 / 9.0;
	image_height = (int)(image_width / aspect_ratio);
	// Viewport
	viewport_height = 2.0;
	viewport_width = aspect_ratio * viewport_height;
	// Camera
	focal_length = 1.0;
	camera_center = (t_vec3){.x = 0, .y = 0, .z = 0};
	// Calculate the vectors across the horizontal and down the vertical viewport edges.
	viewport_u = (t_vec3){.x = viewport_width, .y = 0, .z = 0};
	viewport_v = (t_vec3){.x = 0, .y = viewport_height, .z = 0};
	// Calculate the horizontal and vertical delta vectors from pixel to pixel.
	pixel_delta_u = vec3_scale(viewport_u, 1.0 / (image_width - 1));
	pixel_delta_v = vec3_scale(viewport_v, 1.0 / (image_height - 1));
	// Calculate the location of the upper left pixel.
	viewport_upper_left = vec3_sub(vec3_sub(vec3_sub(camera_center,
				vec3_scale(viewport_u, 0.5)), vec3_scale(viewport_v, 0.5)),
		(t_vec3){.x = 0, .y = 0, .z = focal_length});
	pixel00 = vec3_add(viewport_upper_left, vec3_scale(vec3_add(pixel_delta_u,
				pixel_delta_v), 0.5));
	// Open file
	name = malloc(sizeof(char) * 1000);
	srand(time(NULL));
	sprintf(name, "image_%d.ppm", rand() % 9999);
	file = fopen(name, "w");
	if (file == NULL)
	{
		perror("Error opening file");
		return (1);
	}
	// Write header to the file
	fprintf(file, "P3\n%d %d\n255\n", image_width, image_height);
	// Render
	for (int j = image_height - 1; j >= 0; --j)
	{
		for (int i = 0; i < image_width; ++i)
		{
			pixel_center = vec3_add(vec3_add(pixel00, vec3_scale(pixel_delta_u,
						i)), vec3_scale(pixel_delta_v, j));
			ray_direction = vec3_sub(pixel_center, camera_center);
			ray = (t_ray){.origin = camera_center, .direction = ray_direction};
			color = ray_color(&ray);
			vec3_print(color, file);
		}
	}
	fclose(file);
	return (0);
}
/*
for (int j = 0; j < image_height; ++j)
{
for (int i = 0; i < image_width; ++i)
{
	r = (double)i / (image_width - 1);
	g = (double)j / (image_height - 1);
	b = 0.25;
	color = vec3_add(vec3_add(vec3_scale(red, 1),
			vec3_scale(green, g)),
		vec3_scale(blue, r));
	vec3_print(color, file);
}
}*/