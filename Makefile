raytrace: main.cc camera.h hittable.h hittable_list.h material.h sphere.h rtweekend.h
	g++ -std=c++11 -g -O3 -o raytrace $< -lwebp -I /usr/include/ffmpeg

clean:
	rm raytrace

