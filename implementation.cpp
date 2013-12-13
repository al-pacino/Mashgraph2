#include <stdio.h>
#include "EasyBMP/EasyBMP.h"
#include <vector>
#include <string>
#include <math.h>

using std::vector;
using std::pair;
using std::string;
using std::make_pair;

const double PI = 3.141592653589793;
const float sqrt2 = sqrt(2);

typedef vector<float> TFeature;

template<class Type>
class SimpleMatrix
{
private:
	void freem();
public:
	const int width, height;
	Type **matrix;
	SimpleMatrix(int, int);
	~SimpleMatrix();
	
	void operator=(const SimpleMatrix&);
	SimpleMatrix(const SimpleMatrix &im)
		: matrix(0), width(im.width), height(im.height)
	{
		operator=(im);
	}
};
template<class Type>
SimpleMatrix<Type>::SimpleMatrix(int h, int w)
	: width(w), height(h)
{
	matrix = new Type*[h];
	for(int i = 0; i < height; i++)
		matrix[i] = new Type[width];
}
template<class Type>
SimpleMatrix<Type>::~SimpleMatrix()
{
	freem();
}
template<class Type>
void SimpleMatrix<Type>::operator=(const SimpleMatrix &im)
{
	if(im.matrix == matrix)
		return;
	freem();
	const_cast<int&>(width) = im.width;
	const_cast<int&>(height) = im.height;
	matrix = im.matrix;
	(const_cast<SimpleMatrix<Type>&>(im)).matrix = 0;
}

template<class Type>
inline void SimpleMatrix<Type>::freem()
{
	if(!matrix)
		return;
	for(int i = 0; i < height; i++)
		delete[] matrix[i];
	delete[] matrix;
}

typedef SimpleMatrix<int> IntMatrix;

IntMatrix GetGrayImage(BMP *image)
{
	int w = image->TellWidth(), h = image->TellHeight();
	IntMatrix dest(h, w);

	for(int i = 0; i < h; i++)
	{
		for(int j = 0; j < w; j++)
		{
			RGBApixel pix = image->GetPixel(i, j);
			dest.matrix[i][j] = int(0.299*pix.Red + 0.587*pix.Green + 0.114*pix.Blue);
		}
	}
	return dest;
}

#if 1
TFeature GetFeature(BMP *image)
{
	IntMatrix img = GetGrayImage(image);

	const int default_width = 32;
	const int default_height = 32;
	const int cell_width = 3;
	const int cell_height = 3;
	const int block_width = 2;
	const int block_height = 2;
	const int count_of_segments = 15;

	if(img.width != default_width || img.height != default_height)
		throw string("incorrect width or height. only 32 * 32 bmp availiable");

	int scharr_x[3][3] = {{3, 0, -3},{10, 0, -10},{3, 0, -3}};
	int scharr_y[3][3] = {{3, 10, 3},{0, 0, 0},{-3, -10, -3}};

	int cell_count_x = (default_width-2) / cell_width;
	int cell_count_y = (default_height-2) / cell_height;
	int w = cell_count_x * cell_width + 1, h = cell_count_y * cell_height + 1;
	int size = count_of_segments * cell_count_x * cell_count_y;

	TFeature dest(size, 0.0);

	int dest_shift_y = 0;
	for(int i = 1; i < h; i++)
	{
		int dest_shift_x = 0;
		for(int j = 1; j < w; j++)
		{
			/* calculate gradient x and y */
			int gx = 0, gy = 0;
			for(int ki = -1; ki <= 1; ki++)
			{
				for(int kj = -1; kj <= 1; kj++)
				{
					gx += img.matrix[i+ki][j+kj] * scharr_x[ki+1][kj+1];
					gy += img.matrix[i+ki][j+kj] * scharr_y[ki+1][kj+1];
				}
			}
			/* calculate number of segment (from 0 to count_of_segments-1) */
			int seg = int(floor(count_of_segments*0.5*(1.0+atan2(gy, gx)/PI)));
			if(seg >= count_of_segments)
				seg = count_of_segments-1;
			/* calculate module */
			float mod = sqrt(gx*gx + gy*gy);
			/* hystogram */
			dest[dest_shift_y + dest_shift_x + seg] += mod;
			/* increment shift x */
			if(j % cell_width == 0)
				dest_shift_x += count_of_segments;
		}
		/* increment shift x */
		if(i % cell_height == 0)
			dest_shift_y += count_of_segments * cell_count_x;
	}

	int block_count_x = cell_count_x / block_width;
	int block_count_y = cell_count_y / block_height;
	int sy = 0, sx = 0, step_sx = block_width * count_of_segments;
	int step_ci = count_of_segments * cell_count_x;
	int step_sy = block_height * step_ci;

	for(int i = 0; i < block_count_y; i++, sy += step_sy)
	{
		sx = 0;
		for(int j = 0; j < block_count_x; j++, sx += step_sx)
		{
			double sum = 0.0;
			for(int ci = 0; ci < block_height; ci++)
			{
				for(int cj = 0; cj < step_sx; cj++)
				{
					float v = dest[sy + sx + ci*step_ci + cj];
					sum += v*v;
				}
			}
			if(sum == 0.0)
				continue;
			sum = sqrt(sum);
			for(int ci = 0; ci < block_height; ci++)
			{
				for(int cj = 0; cj < step_sx; cj++)
				{
					dest[sy + sx + ci*step_ci + cj] /= sum;
				}
			}
		}
	}
	return dest;
}
#else
TFeature GetFeature(BMP *image)
{
	IntMatrix img = GetGrayImage(image);

	const int default_width = 32;
	const int default_height = 32;
	const int cell_width = 3;
	const int cell_height = 3;
	const int count_of_segments = 16;

	if(img.width != default_width || img.height != default_height)
		throw string("incorrect width or height. only 32 * 32 bmp availiable");

	int scharr_x[3][3] = {{3, 0, -3},{10, 0, -10},{3, 0, -3}};
	int scharr_y[3][3] = {{3, 10, 3},{0, 0, 0},{-3, -10, -3}};

	int cell_count_x = (default_width-2) / cell_width;
	int cell_count_y = (default_height-2) / cell_height;
	int dsx = (default_width - cell_count_x * cell_width) / 2;
	int dsy = (default_height - cell_count_y * cell_height) / 2;

	int features_length = cell_count_x * cell_count_y * count_of_segments;
	TFeature feature(features_length);

	for(int k = 0, s = 0; k < features_length; k += count_of_segments, s++)
	{
		int sy = dsy + (s / cell_count_x) * cell_height;
		int sx = dsx + (s % cell_count_x) * cell_width;
		for(int i = 0; i < cell_height; i++)
		{
			for(int j = 0; j < cell_width; j++)
			{
				int gr_x = 0, gr_y = 0;
				for(int ki = -1; ki <= 1; ki++)
				{
					for(int kj = -1; kj <= 1; kj++)
					{
						int v = img.matrix[sy+i+ki][sx+j+kj];
						gr_x += v * scharr_x[ki+1][kj+1];
						gr_y += v * scharr_y[ki+1][kj+1];
					}
				}
				float gr = sqrt(gr_y * gr_y + gr_x * gr_x);
				int seg = trunc(count_of_segments*(atan2(gr_y, gr_x)+PI)/(2.0*PI));
				seg = seg >= count_of_segments ? count_of_segments-1 : seg;
				feature[k+seg] = gr;
			}
		}
	}
	TFeature norm_feature(2*(features_length-1));
	int limit = features_length - count_of_segments;
	for(int i = 0; i < limit; i += count_of_segments)
	{
		double sum = 0.0;
		for(int k = 0; k < count_of_segments; k++)
		{
			float v1 = feature[i+k];
			float v2 = feature[i+k+count_of_segments];
			sum += double(v1 + v2);
		}
		if(sum == 0.0)
		{
			for(int k = 0; k < count_of_segments; k++)
			{
				norm_feature.push_back(0);
				norm_feature.push_back(0);
			}
		}
		else
		{
			for(int k = 0; k < count_of_segments; k++)
			{
				norm_feature.push_back(feature[i+k] / sum);
				norm_feature.push_back(feature[i+k+count_of_segments] / sum);
			}
		}
	}
	return norm_feature;
}
#endif

TFeature GetUnlinearFeature(BMP *image)
{
	return GetFeature(image);
	const float L = 0.5f;

	TFeature lf = GetFeature(image);
	TFeature dest(5*lf.size());

	for(TFeature::const_iterator i = lf.begin(); i != lf.end(); ++i)
	{
		float x = *i;
		if(x == 0.0)
		{
			dest.push_back(0.0);
			dest.push_back(0.0);
			dest.push_back(0.0);
			dest.push_back(0.0);
			dest.push_back(0.0);
		}
		else
		{
			float mod = sqrt2 * sqrt(x/(exp(-L*PI)+exp(L*PI)));
			float logxl = L*log(x);
			float re1 = cos(-logxl) * mod;
			float im1 = sin(-logxl) * mod;
			float re3 = cos(L*logxl) * mod;
			float im3 = sin(L*logxl) * mod;
			dest.push_back(re1);
			dest.push_back(im1);
			dest.push_back(mod);
			dest.push_back(re3);
			dest.push_back(im3);
		}
	}
	return dest;
}

