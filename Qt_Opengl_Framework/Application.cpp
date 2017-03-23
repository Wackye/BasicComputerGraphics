#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <set>
#include <stdlib.h>
#include <math.h>
Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene( void )
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	
}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================

void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0 ; i < 3 ; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------
 
///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel should be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = (i*img_width+j)*3; 
			int offset_rgba = (i*img_width+j)*4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k=0; k<3; k++)
				img_data[offset_rgba + k] = gray; //把rgb都丟入gray的顏色, alpha變成白色, 色彩就會變灰階

			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//  
//  Use to find the median, or histogram.   
//   
// 
///////////////////////////////////////////////////////////////////////////////

unsigned char Application::Quick_select(std::vector<unsigned char> gray)
{
	std::sort(gray.begin(), gray.end());
	return(gray[gray.size() / 2 - 1]);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();
	
		for (int i = 0; i < img_height*3; i++)
		{
			for (int j = 0; j < img_width*3; j++)
			{
				img_data[(i * img_width + j)] = (img_data[(i * img_width + j)] / 32) * 32;	
			}
		}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();
	
	
	std::set<int> r, g, b;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = (img_width  * i + j) * 3;
				r.insert(rgb[offset_rgb + rr]);
				r.insert(rgb[offset_rgb + gg]);
				r.insert(rgb[offset_rgb + bb]);
		}
	}



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();
	
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = (i*img_width + j) * 3;
			int offset_rgba = (i*img_width + j) * 4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k = 0; k < 3; k++)
			{ // 直接從gray借來用
				if (gray > 128) img_data[offset_rgba + k] = WHITE;
				else img_data[offset_rgba + k] = BLACK;
			}
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();

	Gray();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = (i*img_width + j) * 4;
			int temp  = (int)(img_data[offset_rgba]) + rand() % 81 - 40;
			if (temp > (int)128)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = WHITE;
			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = BLACK;
			}
			img_data[offset_rgba + 3] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	std::vector<unsigned char> G;
	
	Gray();
	
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			G.push_back(rgb[(i * img_width + j) * 3]);
		}
	}

	unsigned char mid = Quick_select(G);
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = (i * img_width + j) * 4;
			if (img_data[offset_rgb] > mid)
				for(int k = 0; k < 3; k++)	img_data[offset_rgb + k] = WHITE;
			else 
				for(int k = 0; k < 3; k++)	img_data[offset_rgb + k] = BLACK;

			img_data[offset_rgb + 3] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();
	
	Gray();
	double *Pix = new double[img_width * img_height];
	for (int i = 0; i < img_height; i++)
		for (int j = 0; j < img_width; j++)
			Pix[i*img_width + j] = img_data[4 * (i*img_width + j)];

	for (int i = 0; i < img_height-3; i++)
		for (int j = 0; j < img_width-3; j++)
		{
			Pix[i*img_width + j] *= 0.75; //1
			Pix[i*img_width + j+1] *= 0.375; //2
			Pix[i*img_width + j+2] *= 0.625; //3
			Pix[i*img_width + j + 3] *= 0.25; //4
			Pix[(i+1)*img_width + j] *= 0.0625; //5
			Pix[(i+1)*img_width + j+1] *= 1.0; //6
			Pix[(i+1)*img_width + j+2] *= 0.875; //7
			Pix[(i+1)*img_width + j+3] *= 0.4375; //8
			Pix[(i+2)*img_width + j] *= 0.5; //9
			Pix[(i+2)*img_width + j+1] *= 0.8125; //10
			Pix[(i+2)*img_width + j+2] *= 0.9375; //11
			Pix[(i+2)*img_width + j+3] *= 0.125; //12
			Pix[(i+3)*img_width + j] *= 0.1875; //13
			Pix[(i+3)*img_width + j+1] *= 0.5625; //14
			Pix[(i+3)*img_width + j+2] *= 0.3125; //15
			Pix[(i+3)*img_width + j+3] *= 0.6875; //16
		}

	for (int i = 0; i < img_height; i++)
		for (int j = 0; j < img_width; j++)
			for(int k = 0; k < 3; k++)
			{
				Pix[i*img_width + j] >= 32 ? img_data[4 * (i*img_width + j) + k] = WHITE : img_data[4 * (i*img_width + j) + k] = BLACK;
			}

	delete[] Pix;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();
	Gray();

	int *Pix = new int[img_height*img_width]();

	for (int i = 0; i < img_height; i++)
		for (int j = 0; j < img_width; j++)
		{
			int nowIndex = i*img_width + j;
			Pix[nowIndex] = img_data[4 * nowIndex];
		}

	for (int i = 0; i < img_height - 1; i++)
	{
		for (int j = 1; j < img_width - 1 && j > 0; j++)
		{
			int offset = (i*img_width + j);
			int oldPix = Pix[offset];
			int newPix = (oldPix >= 128) ? 255 : 0;
			Pix[offset] = newPix;

			int quant_error = oldPix - newPix;
			/// #region matrix product
			Pix[i    *img_width + j + 1] += double(quant_error * 7 / 16); //right 				
			Pix[(i + 1)*img_width + j - 1] += double(quant_error * 3 / 16); //down-left
			Pix[(i + 1)*img_width + j] += double(quant_error * 5 / 16); //down
			Pix[(i + 1)*img_width + j + 1] += double(quant_error * 1 / 16); //down-right
																			/// #endRegion
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 1; j < img_width; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				img_data[(i*img_width + j) * 4 + k] = Pix[i*img_width + j];
			}
		}
	}



	delete[] Pix;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
	
	int *R = new int[img_height*img_width], *G = new int[img_height*img_width], *B = new int[img_height*img_width];

	for (int i = 0; i < img_height; i++)
		for (int j = 0; j < img_width; j++)
		{
			int nowIndex = i*img_width + j;
			R[nowIndex] = img_data[4 * nowIndex + rr];
			G[nowIndex] = img_data[4 * nowIndex + gg];
			B[nowIndex] = img_data[4 * nowIndex + bb];
		}

	for (int i = 0; i < img_height - 1; i++)
	{
		for (int j = 1; j < img_width - 1 && j > 0; j++)
		{
			int offset = (i*img_width + j);
			int oldR = R[offset + rr];
			int oldG = G[offset + gg];
			int oldB = B[offset + bb];
			int newR = (oldR >= 128) ? 255 : 0;
			int newG = (oldG >= 128) ? 255 : 0;
			int newB = (oldB >= 128) ? 255 : 0;
			R[offset] = newR;
			G[offset] = newG;
			B[offset] = newB;

			int errorR = oldR - newR;
			int errorG = oldG - newG;
			int errorB = oldB - newB;
			/// #region matrix product
			R[i    *img_width + j + 1] += double(errorR * 7 / 16); //right 				
			R[(i + 1)*img_width + j - 1] += double(errorR * 3 / 16); //down-left
			R[(i + 1)*img_width + j] += double(errorR * 5 / 16); //down
			R[(i + 1)*img_width + j + 1] += double(errorR * 1 / 16); //down-right
			
			G[i    *img_width + j + 1] += double(errorG * 7 / 16); //right 				
			G[(i + 1)*img_width + j - 1] += double(errorG * 3 / 16); //down-left
			G[(i + 1)*img_width + j] += double(errorG * 5 / 16); //down
			G[(i + 1)*img_width + j + 1] += double(errorG * 1 / 16); //down-right

			B[i    *img_width + j + 1] += double(errorB * 7 / 16); //right 				
			B[(i + 1)*img_width + j - 1] += double(errorB * 3 / 16); //down-left
			B[(i + 1)*img_width + j] += double(errorB * 5 / 16); //down
			B[(i + 1)*img_width + j + 1] += double(errorB * 1 / 16); //down-right

			/// #endRegion
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 1; j < img_width; j++)
		{
				img_data[(i*img_width + j) * 4 + rr] = R[i*img_width + j];
				img_data[(i*img_width + j) * 4 + gg] = G[i*img_width + j];
				img_data[(i*img_width + j) * 4 + bb] = B[i*img_width + j];
		}
	}



	delete[] R,G,B;

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//-----------------------RGBtoHSV-----------------------
unsigned char* Application::RGB2HSV(unsigned char* RGB)
{
	//H = 0 ~ 360, S = 0 ~ 100, V = 0 ~ 255
	unsigned char *HSV = new unsigned char[img_width * img_height * 3];

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset = 3 * (i*img_width + j);
			unsigned char R = RGB[offset + rr], G = RGB[offset + gg], B = RGB[offset + bb];

			unsigned char rgbMin = R < G ? (R < B ? R : B) : (G < B ? G : B);
			unsigned char rgbMax = R > G ? (R > B ? R : B) : (G > B ? G : B);
			unsigned char H, S, V;
			V = rgbMax;
			if (V == 0)
			{
				HSV[offset] = 0;
				HSV[offset + 1] = 0;
				HSV[offset + 2] = 0;
				continue;
			}

			S = (255 * (rgbMax - rgbMin)) / V;
			if (S == 0)
			{
				HSV[offset] = 0;
				HSV[offset + 1] = 0;
				HSV[offset + 2] = V;
			}
			else
			{
				if (rgbMax == R)
					H =  (43 * (G - B)) / (rgbMax - rgbMin);
				else if (rgbMax == G)
					H = 85 + ( 43 * (B - R)) / (rgbMax - rgbMin);
				else
					H = 171 + ( 43 * (R - G)) / (rgbMax - rgbMin);

					HSV[offset] = H;
					HSV[offset + 1] = S;
					HSV[offset + 2] = V;
			}
		}
	}

	return HSV;
}
	
unsigned char* Application::HSV2RGB(unsigned char* HSV)
{
	unsigned char *RGB = new unsigned char[img_width * img_height * 3];
	for (int i = 0; i < img_height; i++)
	{
		for(int j = 0; j < img_width; j++)
		{
			int offset = (i*img_width + j) * 3;
			unsigned char region, remainder, p, q, t;
			unsigned char H = HSV[offset], S = HSV[offset + 1], V = HSV[offset + 2];
			unsigned R, G, B;
			if (S == 0)
			{
				R = V;
				G = V;
				B = V;
			}

			region = H / 43;
			remainder = (H - (region * 43)) * 6;

			p = (V * (255 - S)) >> 8;
			q = (V * (255 - ((S * remainder) >> 8))) >> 8;
			t = (V * (255 - ((S * (255 - remainder)) >> 8))) >> 8;

			switch (region)
			{
			case 0:
				R = V; G = t; B = p;
				break;
			case 1:
				R = q; G = V; B = p;
				break;
			case 2:
				R = p; G = V; B = t;
				break;
			case 3:
				R = p; G = q; B = V;
				break;
			case 4:
				R = t; G = p; B = V;
				break;
			default:
				R = V; G = p; B = q;
				break;
			}
			RGB[offset + rr] = R, RGB[offset + bb] = B, RGB[offset + gg] = G;
		}
	}
	return RGB;
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering( double filter[][5] )
{
	unsigned char *rgb = this->To_RGB();

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::filtering( double **filter, int n )
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	unsigned char *rgb = this->To_RGB();
	unsigned char *hsv;
	hsv = RGB2HSV(rgb);
	//for(int t = 0; t < 1; t++)
		for (int i = 1; i < img_height - 1; i++)
		{
			for (int j = 1; j < img_width - 1; j++)
			{
				int offset_hsv = 3 * (i*img_width + j);
				// 			if (i-1 >= 0) //如果有上方for (int k = -1; k < 2; k++)
				int sumS = 0, sumV = 0;
				for (int k = -1; k <= 1; k++)
				{
					sumS += hsv[((i - 1)*img_width + j + k) * 3 + 1] + hsv[(i*img_width + j + k) * 3 + 1] + hsv[((i + 1)*img_width + j + k) * 3 + 1];
					//				S    = hsv[上]							+ hsv[中]					 + hsv[下]
					sumV += hsv[((i - 1)*img_width + j + k) * 3 + 2] + hsv[(i*img_width + j + k) * 3 + 2] + hsv[((i + 1)*img_width + j + k) * 3 + 2];
					//				V	 = 
				}
				hsv[offset_hsv + 1] = sumS/9;
				hsv[offset_hsv + 2] = sumV/9;
				// 			{
				// 				tmpS += hsv[3 * ((i - 1)*img_width + j) + 1]; //上
				// 				tmpV += hsv[3 * ((i - 1)*img_width + j) + 2];
				// 				cnt++;
				// 
				// 				if(j-1 >= 0 && j < img_width-1) //如果有左~中間
				// 				{
				// 					tmpS += hsv[3 * ((i-1)*img_width + j-1) + 1]; //左上
				// 					tmpS += hsv[3 * ((i)*img_width + j-1) + 1]; //右上
				// 				}
				// 				else if(j+1 < img_width )
				// 				}
				// 			}
			}
		}
		rgb = HSV2RGB(hsv);
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				
		  		for (int k = 0; k < 3; k++)
				{
					img_data[4 * (img_width*i + j) + k] = rgb[(i*img_width + j)*3 + k ];
				}
			}
		}
		mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
		renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	unsigned char *rgb = this->To_RGB();
	unsigned char *hsv;
	hsv = RGB2HSV(rgb);
	for (int i = 0+2; i < img_height - 2; i++)
		{
			for (int j = 0+2; j < img_width - 2; j++)
			{
				int offset_hsv = 3 * (i*img_width + j);
				// 			if (i-1 >= 0) //如果有上方for (int k = -1; k < 2; k++)
				int sumS = 0, sumV = 0;
				for (int h = -2; h <= 2; h++)
				{
					for (int v = -2; v <= 2; v++)
					{	//horizontal, vertical
						int weightH = Combination(2 + 2, abs(h+2)), weightV = Combination(2 + 2, abs(v+2));
						sumS += weightH * weightV * (int)hsv[((i + h)*img_width + j + v)*3 + 1];
						sumV += weightH * weightV * (int)hsv[((i + h)*img_width + j + v)*3 + 2];
					}
				}
				hsv[offset_hsv + 1] = sumS / 256;
				hsv[offset_hsv + 2] = sumV / 256;
			}
		}
		rgb = HSV2RGB(hsv);
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					img_data[4 * (img_width*i + j) + k] = rgb[(i*img_width + j) * 3 + k];
				}
			}
		}
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();

}
int Application::Combination(int m, int n) //Cm取n
{
	int sumM = 1, sumN = 1, sumMN = 1;
	for (int i = 1; i <= m; i++)
		sumM *= i;
	for (int i = 1; i <= n; i++)
		sumN *= i;
	for (int i = 1; i <= m - n; i++)
		sumMN *= i;

	return sumM / (sumN * sumMN);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N( int N )
{
	if (N % 2 == 0) exit(0);
	unsigned char *rgb = this->To_RGB();
	unsigned char *hsv;
	hsv = RGB2HSV(rgb);
		for (int i = 0 + N; i < img_height - N; i++)
		{
			for (int j = 0 + N; j < img_width - N; j++)
			{
				int offset_hsv = 3 * (i*img_width + j);
				// 			if (i-1 >= 0) //如果有上方for (int k = -1; k < 2; k++)
				int sumS = 0, sumV = 0, sum = 0;
				for (int h = -N/2; h <= N/2; h++)
				{
					for (int v = -N/2; v <= N/2; v++)
					{	//horizontal, vertical
						int weightH = Combination(N-1, h + N/2), weightV = Combination(N-1, v + N / 2);
						sum += weightH * weightV;
						sumS += weightH * weightV * (int)hsv[((i + h)*img_width + j + v) * 3 + 1];
						sumV += weightH * weightV * (int)hsv[((i + h)*img_width + j + v) * 3 + 2];
					}
				}
				hsv[offset_hsv + 1] = sumS / sum;
				hsv[offset_hsv + 2] = sumV / sum;
			}
		}
		rgb = HSV2RGB(hsv);
		for (int i = 0; i < img_height; i++)
		{
			for (int j = 0; j < img_width; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					img_data[4 * (img_width*i + j) + k] = rgb[(i*img_width + j) * 3 + k];
				}
			}
		}
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	 
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////0
void Application::Half_Size()
{
	unsigned char *newRGBA = new unsigned char[img_height*img_width]; //img_height/2*img_width/2 * 4
	size_t count = 0;
	int sum[4] = { 0 };
	for (int i = 0; i < img_height - 1; i+=2)
	{
		for (int j = 0; j < img_width - 1; j+=2)
		{
			for (int k = 0; k < 4; k++)
			{//B G R A //沒有B圖層
				sum[k] = (int)(img_data[(i*img_width + j) * 4 + k] + img_data[(i*img_width + j + 1) * 4 + k] + img_data[((i + 1)*img_width + j) * 4 + k] + img_data[((i + 1)*img_width + j + 1) * 4 + k]) / 4;
			}
			//(i*(img_width)/2 + j/2)*4
				newRGBA[i*img_width + j*2 + 0] = sum[0];
				newRGBA[i*img_width + j*2 + 1] = sum[1];
				newRGBA[i*img_width + j*2 + 2] = sum[2];
				newRGBA[i*img_width + j*2 + 3] = sum[3];
		}
	}
	
	mImageDst = QImage(newRGBA, img_width/2, img_height/2, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}



