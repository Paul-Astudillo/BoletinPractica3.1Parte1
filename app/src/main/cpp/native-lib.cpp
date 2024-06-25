#include <jni.h>
#include <opencv2/opencv.hpp>
#include <android/bitmap.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <android/log.h>

using namespace cv;
using namespace std;

void bitmapToMat(JNIEnv * env, jobject bitmap, cv::Mat &dst, jboolean needUnPremultiplyAlpha) {
    AndroidBitmapInfo info;
    void* pixels = 0;
    try {
        CV_Assert( AndroidBitmap_getInfo(env, bitmap, &info) >= 0 );
        CV_Assert( info.format == ANDROID_BITMAP_FORMAT_RGBA_8888 ||
                   info.format == ANDROID_BITMAP_FORMAT_RGB_565 );
        CV_Assert( AndroidBitmap_lockPixels(env, bitmap, &pixels) >= 0 );
        CV_Assert( pixels );
        dst.create(info.height, info.width, CV_8UC4);
        if( info.format == ANDROID_BITMAP_FORMAT_RGBA_8888 ) {
            cv::Mat tmp(info.height, info.width, CV_8UC4, pixels);
            if(needUnPremultiplyAlpha) cvtColor(tmp, dst, cv::COLOR_mRGBA2RGBA);
            else tmp.copyTo(dst);
        } else {
            // info.format == ANDROID_BITMAP_FORMAT_RGB_565
            cv::Mat tmp(info.height, info.width, CV_8UC2, pixels);
            cvtColor(tmp, dst, cv::COLOR_BGR5652RGBA);
        }
        AndroidBitmap_unlockPixels(env, bitmap);
        return;
    } catch(const cv::Exception& e) {
        AndroidBitmap_unlockPixels(env, bitmap);
        jclass je = env->FindClass("java/lang/Exception");
        env->ThrowNew(je, e.what());
        return;
    } catch (...) {
        AndroidBitmap_unlockPixels(env, bitmap);
        jclass je = env->FindClass("java/lang/Exception");
        env->ThrowNew(je, "Unknown exception in JNI code {nBitmapToMat}");
        return;
    }
}
// Momentos Hu de referencia
vector<double> MOMENTOS_HU_CIRCULO = {6.40926530e-04, 1.86979797e-08, 2.00214655e-12, 2.26706182e-14, 1.34407811e-27, 1.43894564e-19, -3.07094768e-28};
vector<double> MOMENTOS_HU_TRIANGULO = {8.06510150e-04, 5.73788071e-08, 3.48045489e-10, 1.10871558e-11, 8.84029050e-22, 2.47656459e-15, 4.38685134e-22};
vector<double> MOMENTOS_HU_CUADRADO = {6.85473402e-04, 3.44629921e-08, 4.16311748e-12, 1.13729153e-13, 3.79966091e-26, -1.50169729e-18, -7.60574069e-26};

// Función para obtener los momentos Hu de una imagen
vector<double> get_hu_moments(const Mat& src) {
    // Convertir la imagen a HSV
    Mat hsv;
    cvtColor(src, hsv, COLOR_BGR2HSV);

    // Rango para binarizar la imagen (estos valores pueden ajustarse según sea necesario)
    Scalar lower_bound(107, 20, 20);
    Scalar upper_bound(140, 255, 255);

    // Binarizar la imagen
    Mat binary;
    inRange(hsv, lower_bound, upper_bound, binary);

    // Aplicar filtro de mediana
    medianBlur(binary, binary, 3);

    // Convertir la imagen a color para poder usar floodFill
    Mat colored;
    cvtColor(binary, colored, COLOR_GRAY2BGR);

    // Punto inicial
    Point seed_point(colored.cols / 2, colored.rows / 2);

    // Color de relleno blanco
    Scalar color_relleno(255, 255, 255);

    // Aplicar floodFill
    floodFill(colored, seed_point, color_relleno);

    // Convertir la imagen de nuevo a escala de grises
    Mat relleno;
    cvtColor(colored, relleno, COLOR_BGR2GRAY);

    // Calcular los momentos de la imagen
    Moments moments = cv::moments(relleno);
    double hu[7];
    HuMoments(moments, hu);

    // Convertir a vector de double
    vector<double> hu_moments(hu, hu + 7);

    // Imprimir los momentos Hu para depuración
    __android_log_print(ANDROID_LOG_DEBUG, "MOMENTOS_HU", "Hu Moments: %f, %f, %f, %f, %f, %f, %f", hu[0], hu[1], hu[2], hu[3], hu[4], hu[5], hu[6]);

    return hu_moments;
}

// Función para calcular la distancia euclidiana entre dos vectores
double euclidean_distance(const vector<double>& vec1, const vector<double>& vec2) {
    double sum = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        sum += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(sum);
}

// Función para clasificar la figura basada en los momentos Hu
string classify_shape(const Mat& src, const vector<double>& square_mean_hu, const vector<double>& triangle_mean_hu, const vector<double>& circle_mean_hu) {
    vector<double> hu_moments = get_hu_moments(src);

    map<string, double> distances;
    distances["cuadrado"] = euclidean_distance(hu_moments, square_mean_hu);
    distances["triangulo"] = euclidean_distance(hu_moments, triangle_mean_hu);
    distances["circulo"] = euclidean_distance(hu_moments, circle_mean_hu);

    string predicted_shape;
    double min_distance = numeric_limits<double>::max();
    for (const auto& pair : distances) {
        if (pair.second < min_distance) {
            min_distance = pair.second;
            predicted_shape = pair.first;
        }
    }

    // Imprimir las distancias para depuración
    __android_log_print(ANDROID_LOG_DEBUG, "CLASSIFICATION", "Distancias: Cuadrado: %f, Triángulo: %f, Círculo: %f", distances["cuadrado"], distances["triangulo"], distances["circulo"]);

    return predicted_shape;
}

// Función JNIEXPORT para clasificar una figura desde una imagen Bitmap
extern "C" JNIEXPORT jstring JNICALL
Java_com_example_boletinpractica3_MainActivity_calcularHuMoments
        (JNIEnv* env, jobject obj, jobject fotoObj) {
    Mat fotoMat;
    bitmapToMat(env, fotoObj, fotoMat, false);

    if (fotoMat.empty()) {
        return env->NewStringUTF("Error loading image");
    }

    // Mejorar la calidad de la imagen borrosa usando filtro de desenfoque gaussiano inverso
    Mat blurred, highContrast;
    GaussianBlur(fotoMat, blurred, Size(0, 0), 3);
    addWeighted(fotoMat, 1.5, blurred, -0.5, 0, highContrast);

    try {
        string predicted_shape = classify_shape(highContrast, MOMENTOS_HU_CUADRADO, MOMENTOS_HU_TRIANGULO, MOMENTOS_HU_CIRCULO);
        return env->NewStringUTF(predicted_shape.c_str());
    } catch (const exception& e) {
        return env->NewStringUTF("Error classifying shape");
    }
}


//momneto zernike
#define PI 3.14159265358979323846264338328
#define MAX_L 32
#define MAX_D 15
#define MAX_Z 72
#define MAX_LUT 240

void mb_Znl(double *X, double *Y, double *P, int size, double D, double m10_m00, double m01_m00, double R, double psum, double *zvalues, long *output_size) {
    static double LUT[MAX_LUT];
    static int n_s[MAX_Z], l_s[MAX_Z];
    static char init_lut = 0;

    double x, y, p;
    int i, m, theZ, theLUT, numZ = 0;
    int n = 0, l = 0;
    complex<double> sum[MAX_Z];
    complex<double> Vnl;

    assert(D == MAX_D);

    if (!init_lut) {
        theZ = 0;
        theLUT = 0;
        for (n = 0; n <= MAX_D; n++) {
            for (l = 0; l <= n; l++) {
                if ((n - l) % 2 == 0) {
                    for (m = 0; m <= (n - l) / 2; m++) {
                        LUT[theLUT] = pow((double)-1.0, (double)m) * ((long double)tgamma(n - m + 1) / ((long double)tgamma(m + 1) * (long double)tgamma((n - 2 * m + l) / 2 + 1) * (long double)tgamma((n - 2 * m - l) / 2 + 1)));
                        theLUT++;
                    }
                    n_s[theZ] = n;
                    l_s[theZ] = l;
                    theZ++;
                }
            }
        }
        init_lut = 1;
    }

    for (n = 0; n <= D; n++) {
        for (l = 0; l <= n; l++) {
            if ((n - l) % 2 == 0) {
                sum[numZ] = complex<double>(0.0, 0.0);
                numZ++;
            }
        }
    }

    for (i = 0; i < size; i++) {
        x = (X[i] - m10_m00) / R;
        y = (Y[i] - m01_m00) / R;
        double sqr_x2y2 = sqrt(x * x + y * y);
        if (sqr_x2y2 > 1.0) continue;

        p = P[i] / psum;
        double atan2yx = atan2(y, x);
        theLUT = 0;
        for (theZ = 0; theZ < numZ; theZ++) {
            n = n_s[theZ];
            l = l_s[theZ];
            Vnl = complex<double>(0.0, 0.0);
            for (m = 0; m <= (n - l) / 2; m++) {
                Vnl += (polar(1.0, l * atan2yx) * LUT[theLUT] * pow(sqr_x2y2, (double)(n - 2 * m)));
                theLUT++;
            }
            sum[theZ] += (conj(Vnl) * p);
        }
    }

    for (theZ = 0; theZ < numZ; theZ++) {
        sum[theZ] *= ((n_s[theZ] + 1) / PI);
        double preal = real(sum[theZ]);
        double pimag = imag(sum[theZ]);
        zvalues[theZ] = fabs(sqrt(preal * preal + pimag * pimag));
    }

    *output_size = numZ;
}

void mb_zernike2D(const Mat &Im, double D, double R, double *zvalues, long *output_size) {
    int rows = Im.rows, cols = Im.cols;
    if (D <= 0) D = 15;
    if (R <= 0) R = rows / 2.0;

    double *Y = new double[rows * cols];
    double *X = new double[rows * cols];
    double *P = new double[rows * cols];
    double psum = 0.0;
    double moment10 = 0.0, moment00 = 0.0, moment01 = 0.0;

    int size = 0;
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            double intensity = Im.at<uchar>(y, x);
            if (intensity != 0) {
                Y[size] = y + 1;
                X[size] = x + 1;
                P[size] = intensity;
                psum += intensity;
                size++;
            }
            moment10 += (x + 1) * intensity;
            moment00 += intensity;
            moment01 += (y + 1) * intensity;
        }
    }

    double m10_m00 = moment10 / moment00;
    double m01_m00 = moment01 / moment00;
    mb_Znl(X, Y, P, size, D, m10_m00, m01_m00, R, psum, zvalues, output_size);

    delete[] Y;
    delete[] X;
    delete[] P;
}

void calculateZernikeMoments(const Mat& image, double* zvalues, long* output_size) {
    double D = 15;
    double R = image.rows / 2.0;
    mb_zernike2D(image, D, R, zvalues, output_size);
}

double distanciaEuclidea(double* momentosZernike, double* referencia, int size) {
    double suma = 0;
    for (int i = 0; i < size; ++i) {
        suma += pow(referencia[i] - momentosZernike[i], 2);
    }
    return sqrt(suma);
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_example_boletinpractica3_MainActivity_calcularZernikeMoments
        (JNIEnv* env, jobject obj, jobject fotoObj) {


    double triangulo_zernike[MAX_Z] = {0.30313, 0.00707924, 0.355111, 0.059828, 0.0353174, 0.10038, 0.147648, 0.0684501, 0.0243417, 0.0296546, 0.121543, 0.0461651, 0.0726276, 0.0469767, 0.0343686, 0.058603, 0.0311033, 0.0666696, 0.0555346, 0.0269396, 0.0711033, 0.0347564, 0.0328548, 0.0570468, 0.0389602, 0.0350478, 0.0450958, 0.0424141, 0.0319205, 0.0359617, 0.0469623, 0.0411652, 0.0279026, 0.0347637, 0.0426209, 0.0283167, 0.0298993, 0.0409201, 0.0330015, 0.0258557, 0.032692, 0.0287651, 0.0401726, 0.0375587, 0.0241437, 0.0349209, 0.0388579, 0.0316572, 0.0293563, 0.0277678, 0.0342242, 0.0329382, 0.0264304, 0.0248166, 0.0295647, 0.023145, 0.0335905, 0.0327891, 0.0263654, 0.0282508, 0.0340594, 0.0268838, 0.0266347, 0.0252813, 0.026012, 0.0318314, 0.0379678, 0.0227345, 0.0269958, 0.0271511, 0.0243807, 0.0222122};
    double cuadrado_zernike[MAX_Z] = {0.29503, 0.00282691, 0.209036, 0.0671205, 0.0094695, 0.0176218, 0.0679464, 0.0589661, 0.059665, 0.00918045, 0.0108105, 0.0140892, 0.0488247, 0.0418737, 0.0404804, 0.0362726, 0.0104449, 0.00963323, 0.0116442, 0.0145017, 0.0333797, 0.0387007, 0.033995, 0.0330836, 0.0266531, 0.00970724, 0.0101832, 0.0108972, 0.0125275, 0.0139, 0.0303647, 0.0277955, 0.0291986, 0.025562, 0.0230942, 0.0202858, 0.00990155, 0.00917472, 0.0111056, 0.010475, 0.0110074, 0.0121496, 0.0220006, 0.0240917, 0.0228216, 0.0240465, 0.0210916, 0.0200525, 0.0170988, 0.00918529, 0.0093346, 0.010071, 0.0109331, 0.0112188, 0.0114988, 0.0120059, 0.0213247, 0.0207907, 0.0205473, 0.0203458, 0.0199367, 0.0175176, 0.016988, 0.014419, 0.00987116, 0.00975148, 0.00906047, 0.0109068, 0.00991287, 0.0107677, 0.0099922, 0.0109143};
    double circulo_zernike[MAX_Z] = {0.306341, 0.000439282, 0.3156, 0.0643414, 0.0101655, 0.0142876, 0.108442, 0.0670122, 0.0241566, 0.0111446, 0.0168795, 0.00865969, 0.0960931, 0.0512503, 0.0246005, 0.0127276, 0.0165231, 0.0153324, 0.0124217, 0.00577756, 0.0745371, 0.0545346, 0.0282807, 0.0145599, 0.00852356, 0.0193681, 0.0153204, 0.0129941, 0.00844591, 0.00444148, 0.0563957, 0.0512041, 0.030028, 0.0178813, 0.0109007, 0.00650399, 0.0184655, 0.0171394, 0.0145217, 0.0104351, 0.00672881, 0.00321962, 0.0609544, 0.0461459, 0.0326647, 0.0185547, 0.011622, 0.00879158, 0.00494632, 0.0222841, 0.0171215, 0.0151037, 0.0117948, 0.00803939, 0.00548741, 0.00266032, 0.0422501, 0.0464891, 0.0328822, 0.0222146, 0.01368, 0.00949507, 0.00701259, 0.00413907, 0.0177726, 0.0184439, 0.0150959, 0.0117539, 0.00855819, 0.0062942, 0.00490648, 0.00204912};


    Mat fotoMat;
    bitmapToMat(env, fotoObj, fotoMat, false);

    // Convertir la imagen a HSV
    Mat hsv;
    cvtColor(fotoMat, hsv, COLOR_BGR2HSV);

    // Definir límites de binarización en HSV
    Scalar lowerBound(107, 20, 20); // límite inferior
    Scalar upperBound(140, 255, 255); // límite superior

    // Binarizar la imagen en HSV
    Mat binary;
    inRange(hsv, lowerBound, upperBound, binary);

    // Aplicar un filtro de mediana para suavizar la imagen binarizada
    medianBlur(binary, binary, 3);

    // Convertir la imagen binarizada a color para poder usar floodFill
    Mat colored;
    cvtColor(binary, colored, COLOR_GRAY2BGR);

    // Punto inicial para floodFill
    Point seedPoint(colored.cols / 2, colored.rows / 2);

    // Color de relleno blanco
    Scalar colorRelleno(255, 255, 255);

    // Aplicar flood fill
    floodFill(colored, seedPoint, colorRelleno);

    // Convertir la imagen rellenada a escala de grises
    Mat relleno;
    cvtColor(colored, relleno, COLOR_BGR2GRAY);

    // Calcular momentos de Zernike
    double zernikeMoments[MAX_Z];
    long output_size;
    calculateZernikeMoments(relleno, zernikeMoments, &output_size);

    // Calcular la distancia euclidiana a cada conjunto de momentos promedio
    double distTriangulo = distanciaEuclidea(zernikeMoments, triangulo_zernike, output_size);
    double distCuadrado = distanciaEuclidea(zernikeMoments, cuadrado_zernike, output_size);
    double distCirculo = distanciaEuclidea(zernikeMoments, circulo_zernike, output_size);

    // Determinar la figura con la menor distancia
    string resultado;
    if (distTriangulo < distCuadrado && distTriangulo < distCirculo) {
        resultado = "triangulo";
    } else if (distCuadrado < distTriangulo && distCuadrado < distCirculo) {
        resultado = "cuadrado";
    } else {
        resultado = "círculo";
    }

    return env->NewStringUTF(resultado.c_str());

}