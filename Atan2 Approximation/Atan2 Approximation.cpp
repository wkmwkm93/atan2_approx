// Atan2 Approximation.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>	
using namespace std;

#define M_PI       3.14159265358979323846   // pi
#define M_PI_2     1.57079632679489661923   // pi/2
#define	SQRT_3	1.7320508075688772935274463415059  //1.7320508

#define A1  0.99997726
#define A3  -0.33262347
#define A5  0.19354346
#define A7  -0.11643287
#define A9  0.05265332
#define A11 -0.01172120

void atan2_baseline(size_t num_points, const float* ys, const float* xs, float* out);

void randomFloat(float &dY, float &dX);
inline float atan_scalar_approximation(float x);
inline float atan_scalar_approximation_global(float x);
void atan2_auto_1(size_t num_points, const float* ys, const float* xs, float* out);
void atan2_auto_2(size_t num_points, const float* ys, const float* xs, float* out);
void atan2_auto_2_onebyone(const float ys, const float xs, float out);
double GetATANDeg(double dY, double dX);

const float ATAN_LOOKUP[899] =
{
	0.00175f,0.00349f,0.00524f,0.00698f,0.00873f,0.01047f,0.01222f,0.01396f,0.01571f,0.01746f,
	0.01920f,0.02095f,0.02269f,0.02444f,0.02619f,0.02793f,0.02968f,0.03143f,0.03317f,0.03492f,
	0.03667f,0.03842f,0.04016f,0.04191f,0.04366f,0.04541f,0.04716f,0.04891f,0.05066f,0.05241f,
	0.05416f,0.05591f,0.05766f,0.05941f,0.06116f,0.06291f,0.06467f,0.06642f,0.06817f,0.06993f,
	0.07168f,0.07344f,0.07519f,0.07695f,0.07870f,0.08046f,0.08221f,0.08397f,0.08573f,0.08749f,
	0.08925f,0.09101f,0.09277f,0.09453f,0.09629f,0.09805f,0.09981f,0.10158f,0.10334f,0.10510f,
	0.10687f,0.10863f,0.11040f,0.11217f,0.11394f,0.11570f,0.11747f,0.11924f,0.12101f,0.12278f,
	0.12456f,0.12633f,0.12810f,0.12988f,0.13165f,0.13343f,0.13521f,0.13698f,0.13876f,0.14054f,
	0.14232f,0.14410f,0.14588f,0.14767f,0.14945f,0.15124f,0.15302f,0.15481f,0.15660f,0.15838f,
	0.16017f,0.16196f,0.16376f,0.16555f,0.16734f,0.16914f,0.17093f,0.17273f,0.17453f,0.17633f,
	0.17813f,0.17993f,0.18173f,0.18353f,0.18534f,0.18714f,0.18895f,0.19076f,0.19257f,0.19438f,
	0.19619f,0.19801f,0.19982f,0.20164f,0.20345f,0.20527f,0.20709f,0.20891f,0.21073f,0.21256f,
	0.21438f,0.21621f,0.21804f,0.21986f,0.22169f,0.22353f,0.22536f,0.22719f,0.22903f,0.23087f,
	0.23271f,0.23455f,0.23639f,0.23823f,0.24008f,0.24193f,0.24377f,0.24562f,0.24747f,0.24933f,
	0.25118f,0.25304f,0.25490f,0.25676f,0.25862f,0.26048f,0.26235f,0.26421f,0.26608f,0.26795f,
	0.26982f,0.27169f,0.27357f,0.27545f,0.27732f,0.27921f,0.28109f,0.28297f,0.28486f,0.28675f,
	0.28864f,0.29053f,0.29242f,0.29432f,0.29621f,0.29811f,0.30001f,0.30192f,0.30382f,0.30573f,
	0.30764f,0.30955f,0.31147f,0.31338f,0.31530f,0.31722f,0.31914f,0.32106f,0.32299f,0.32492f,
	0.32685f,0.32878f,0.33072f,0.33266f,0.33460f,0.33654f,0.33848f,0.34043f,0.34238f,0.34433f,
	0.34628f,0.34824f,0.35020f,0.35216f,0.35412f,0.35608f,0.35805f,0.36002f,0.36199f,0.36397f,
	0.36595f,0.36793f,0.36991f,0.37190f,0.37388f,0.37588f,0.37787f,0.37986f,0.38186f,0.38386f,
	0.38587f,0.38787f,0.38988f,0.39190f,0.39391f,0.39593f,0.39795f,0.39997f,0.40200f,0.40403f,
	0.40606f,0.40809f,0.41013f,0.41217f,0.41421f,0.41626f,0.41831f,0.42036f,0.42242f,0.42447f,
	0.42654f,0.42860f,0.43067f,0.43274f,0.43481f,0.43689f,0.43897f,0.44105f,0.44314f,0.44523f,
	0.44732f,0.44942f,0.45152f,0.45362f,0.45573f,0.45784f,0.45995f,0.46206f,0.46418f,0.46631f,
	0.46843f,0.47056f,0.47270f,0.47483f,0.47698f,0.47912f,0.48127f,0.48342f,0.48557f,0.48773f,
	0.48989f,0.49206f,0.49423f,0.49640f,0.49858f,0.50076f,0.50295f,0.50514f,0.50733f,0.50953f,
	0.51173f,0.51393f,0.51614f,0.51835f,0.52057f,0.52279f,0.52501f,0.52724f,0.52947f,0.53171f,
	0.53395f,0.53620f,0.53844f,0.54070f,0.54296f,0.54522f,0.54748f,0.54975f,0.55203f,0.55431f,
	0.55659f,0.55888f,0.56117f,0.56347f,0.56577f,0.56808f,0.57039f,0.57271f,0.57503f,0.57735f,
	0.57968f,0.58201f,0.58435f,0.58670f,0.58905f,0.59140f,0.59376f,0.59612f,0.59849f,0.60086f,
	0.60324f,0.60562f,0.60801f,0.61040f,0.61280f,0.61520f,0.61761f,0.62003f,0.62245f,0.62487f,
	0.62730f,0.62973f,0.63217f,0.63462f,0.63707f,0.63953f,0.64199f,0.64446f,0.64693f,0.64941f,
	0.65189f,0.65438f,0.65688f,0.65938f,0.66189f,0.66440f,0.66692f,0.66944f,0.67197f,0.67451f,
	0.67705f,0.67960f,0.68215f,0.68471f,0.68728f,0.68985f,0.69243f,0.69502f,0.69761f,0.70021f,
	0.70281f,0.70542f,0.70804f,0.71066f,0.71329f,0.71593f,0.71857f,0.72122f,0.72388f,0.72654f,
	0.72921f,0.73189f,0.73457f,0.73726f,0.73996f,0.74267f,0.74538f,0.74810f,0.75082f,0.75355f,
	0.75629f,0.75904f,0.76180f,0.76456f,0.76733f,0.77010f,0.77289f,0.77568f,0.77848f,0.78129f,
	0.78410f,0.78692f,0.78975f,0.79259f,0.79544f,0.79829f,0.80115f,0.80402f,0.80690f,0.80978f,
	0.81268f,0.81558f,0.81849f,0.82141f,0.82434f,0.82727f,0.83022f,0.83317f,0.83613f,0.83910f,
	0.84208f,0.84507f,0.84806f,0.85107f,0.85408f,0.85710f,0.86014f,0.86318f,0.86623f,0.86929f,
	0.87236f,0.87543f,0.87852f,0.88162f,0.88473f,0.88784f,0.89097f,0.89410f,0.89725f,0.90040f,
	0.90357f,0.90674f,0.90993f,0.91313f,0.91633f,0.91955f,0.92277f,0.92601f,0.92926f,0.93252f,
	0.93578f,0.93906f,0.94235f,0.94565f,0.94896f,0.95229f,0.95562f,0.95897f,0.96232f,0.96569f,
	0.96907f,0.97246f,0.97586f,0.97927f,0.98270f,0.98613f,0.98958f,0.99304f,0.99652f,1.00000f,
	1.00350f,1.00701f,1.01053f,1.01406f,1.01761f,1.02117f,1.02474f,1.02832f,1.03192f,1.03553f,
	1.03915f,1.04279f,1.04644f,1.05010f,1.05378f,1.05747f,1.06117f,1.06489f,1.06862f,1.07237f,
	1.07613f,1.07990f,1.08369f,1.08749f,1.09131f,1.09514f,1.09899f,1.10285f,1.10672f,1.11061f,
	1.11452f,1.11844f,1.12238f,1.12633f,1.13029f,1.13428f,1.13828f,1.14229f,1.14632f,1.15037f,
	1.15443f,1.15851f,1.16261f,1.16672f,1.17085f,1.17500f,1.17916f,1.18334f,1.18754f,1.19175f,
	1.19599f,1.20024f,1.20451f,1.20879f,1.21310f,1.21742f,1.22176f,1.22612f,1.23050f,1.23490f,
	1.23931f,1.24375f,1.24820f,1.25268f,1.25717f,1.26169f,1.26622f,1.27077f,1.27535f,1.27994f,
	1.28456f,1.28919f,1.29385f,1.29853f,1.30323f,1.30795f,1.31269f,1.31745f,1.32224f,1.32704f,
	1.33187f,1.33673f,1.34160f,1.34650f,1.35142f,1.35637f,1.36134f,1.36633f,1.37134f,1.37638f,
	1.38145f,1.38653f,1.39165f,1.39679f,1.40195f,1.40714f,1.41235f,1.41759f,1.42286f,1.42815f,
	1.43347f,1.43881f,1.44418f,1.44958f,1.45501f,1.46046f,1.46595f,1.47146f,1.47699f,1.48256f,
	1.48816f,1.49378f,1.49944f,1.50512f,1.51084f,1.51658f,1.52235f,1.52816f,1.53400f,1.53986f,
	1.54576f,1.55170f,1.55766f,1.56366f,1.56969f,1.57575f,1.58184f,1.58797f,1.59414f,1.60033f,
	1.60657f,1.61283f,1.61914f,1.62548f,1.63185f,1.63826f,1.64471f,1.65120f,1.65772f,1.66428f,
	1.67088f,1.67752f,1.68419f,1.69091f,1.69766f,1.70446f,1.71129f,1.71817f,1.72509f,1.73205f,
	1.73905f,1.74610f,1.75319f,1.76032f,1.76749f,1.77471f,1.78198f,1.78929f,1.79665f,1.80405f,
	1.81150f,1.81899f,1.82654f,1.83413f,1.84177f,1.84946f,1.85720f,1.86499f,1.87283f,1.88073f,
	1.88867f,1.89667f,1.90472f,1.91282f,1.92098f,1.92920f,1.93746f,1.94579f,1.95417f,1.96261f,
	1.97111f,1.97966f,1.98828f,1.99695f,2.00569f,2.01449f,2.02335f,2.03227f,2.04125f,2.05030f,
	2.05942f,2.06860f,2.07785f,2.08716f,2.09654f,2.10600f,2.11552f,2.12511f,2.13477f,2.14451f,
	2.15432f,2.16420f,2.17416f,2.18419f,2.19430f,2.20449f,2.21475f,2.22510f,2.23553f,2.24604f,
	2.25663f,2.26730f,2.27806f,2.28891f,2.29984f,2.31086f,2.32197f,2.33317f,2.34447f,2.35585f,
	2.36733f,2.37891f,2.39058f,2.40235f,2.41421f,2.42618f,2.43825f,2.45043f,2.46270f,2.47509f,
	2.48758f,2.50018f,2.51289f,2.52571f,2.53865f,2.55170f,2.56487f,2.57815f,2.59156f,2.60509f,
	2.61874f,2.63252f,2.64642f,2.66046f,2.67462f,2.68892f,2.70335f,2.71792f,2.73263f,2.74748f,
	2.76247f,2.77761f,2.79289f,2.80833f,2.82391f,2.83965f,2.85555f,2.87161f,2.88783f,2.90421f,
	2.92076f,2.93748f,2.95437f,2.97144f,2.98868f,3.00611f,3.02372f,3.04152f,3.05950f,3.07768f,
	3.09606f,3.11464f,3.13341f,3.15240f,3.17159f,3.19100f,3.21063f,3.23048f,3.25055f,3.27085f,
	3.29139f,3.31216f,3.33317f,3.35443f,3.37594f,3.39771f,3.41973f,3.44202f,3.46458f,3.48741f,
	3.51053f,3.53393f,3.55761f,3.58160f,3.60588f,3.63048f,3.65538f,3.68061f,3.70616f,3.73205f,
	3.75828f,3.78485f,3.81177f,3.83906f,3.86671f,3.89474f,3.92316f,3.95196f,3.98117f,4.01078f,
	4.04081f,4.07127f,4.10216f,4.13350f,4.16530f,4.19756f,4.23030f,4.26352f,4.29724f,4.33148f,
	4.36623f,4.40152f,4.43735f,4.47374f,4.51071f,4.54826f,4.58641f,4.62518f,4.66458f,4.70463f,
	4.74534f,4.78673f,4.82882f,4.87162f,4.91516f,4.95945f,5.00451f,5.05037f,5.09704f,5.14455f,
	5.19293f,5.24218f,5.29235f,5.34345f,5.39552f,5.44857f,5.50264f,5.55777f,5.61397f,5.67128f,
	5.72974f,5.78938f,5.85024f,5.91236f,5.97576f,6.04051f,6.10664f,6.17419f,6.24321f,6.31375f,
	6.38587f,6.45961f,6.53503f,6.61219f,6.69116f,6.77199f,6.85475f,6.93952f,7.02637f,7.11537f,
	7.20661f,7.30018f,7.39616f,7.49465f,7.59575f,7.69957f,7.80622f,7.91582f,8.02848f,8.14435f,
	8.26355f,8.38625f,8.51259f,8.64275f,8.77689f,8.91520f,9.05789f,9.20516f,9.35724f,9.51436f,
	9.67680f,9.84482f,10.01871f,10.19879f,10.38540f,10.57889f,10.77967f,10.98815f,11.20478f,11.43005f,
	11.66450f,11.90868f,12.16324f,12.42883f,12.70620f,12.99616f,13.29957f,13.61741f,13.95072f,14.30067f,
	14.66853f,15.05572f,15.46381f,15.89454f,16.34986f,16.83191f,17.34315f,17.88631f,18.46447f,19.08114f,
	19.74029f,20.44649f,21.20495f,22.02171f,22.90377f,23.85928f,24.89783f,26.03074f,27.27149f,28.63625f,
	30.14462f,31.82052f,33.69351f,35.80055f,38.18846f,40.91741f,44.06611f,47.73950f,52.08067f,57.28996f,
	63.65674f,71.61507f,81.84704f,95.48948f,114.58865f,143.23712f,190.98419f,286.47773f,572.95721f
};

int main()
{
	const int SIZE = 1000;
	const int iter = 10;
	float *x = new float[SIZE];
	float *y = new float[SIZE];
	double *x_d = new double[SIZE];
	double *y_d = new double[SIZE];

	float *res = new float[SIZE];

	const int METHODS = 3;

	for (int i = 0; i < SIZE; i++)
	{
		randomFloat(x[i], y[i]);
		x_d[i] = double(x[i]);
		y_d[i] = double(y[i]);
		res[i] = 0.0;
	}

	/*x[0] = 0; y[0] = 0;
	x[1] = 0; y[1] = 0;
	x[2] = 0; y[2] = 0;
	x[3] = 0; y[3] = 0;

	x_d[0] = 0.0; y_d[0] = 0.0;
	x_d[1] = 0.0; y_d[1] = 0.0;
	x_d[2] = 0.0; y_d[2] = 0.0;
	x_d[3] = 0.0; y_d[3] = 0.0;
*/
	/*y[0] = 500; x[0] = 20;
	y[1] = -500; x[1] = 20;
	y[2] = 500; x[2] = -20;
	y[3] = -500; x[3] = -20;
*/
	chrono::duration<double, std::milli> TotalSpan1;
	chrono::duration<double, std::milli> TotalSpan2;
	chrono::duration<double, std::milli> TotalSpan3;

	for (int j = 0; j < iter; j++)
	{
		auto begin = std::chrono::high_resolution_clock::now();
		atan2_baseline(SIZE, y, x, res);
		auto end = std::chrono::high_resolution_clock::now();
		chrono::duration<double, std::milli> elapsed = end - begin;

		//printf("ATAN2 baseline Iter: %d Total Cycle : %.9f ms.\n", j + 1, elapsed.count());
		TotalSpan1 += elapsed;

		begin = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < SIZE; i++)
		{
			double res = GetATANDeg(y_d[i], x_d[i]);
			if(i == SIZE -1)
				printf("GetATANDeg %.2f / %.2f = %.2f\n", y_d[i], x_d[i], res);
		}

		end = std::chrono::high_resolution_clock::now();
		elapsed = end - begin;
		TotalSpan2 += elapsed;

		//printf("ATAN2 AUTO 1 Iter: %d Total Cycle : %.9f ms.\n", j + 1, elapsed.count());

		begin = std::chrono::high_resolution_clock::now();

		atan2_auto_2(SIZE, y, x, res);

		end = std::chrono::high_resolution_clock::now();
		elapsed = end - begin;
		TotalSpan3 += elapsed;

		//printf("ATAN2 AUTO 2 Iter: %d Total Cycle : %.9f ms.\n", j + 1, elapsed.count());
	}

	//for (int i = 0; i<METHODS; i++)
	printf("Method %d ,Each Func Time: %.9f ns.\n", 1, TotalSpan1 *1000000 / (SIZE*iter));
	printf("Method %d ,Each Func Time: %.9f ns.\n", 2, TotalSpan2 *1000000 / (SIZE*iter));
	printf("Method %d ,Each Func Time: %.9f ns.\n", 3, TotalSpan3 *1000000 / (SIZE*iter));

	getchar();
	return 0;
}

void randomFloat(float &dY, float &dX)
{
	float f1 = rand() % 255;
	float f2 = rand() % 255;
	float f3 = rand() % 255;
	float f4 = rand() % 255;
	float f5 = rand() % 255;
	float f6 = rand() % 255;

	dY = SQRT_3*(f2+ f3 - f5 - f6);
	dX = 2.0*f1 + f2 - f3 - 2.0*f4 - f5 + f6;
}

void atan2_baseline(size_t num_points, const float* ys, const float* xs, float* out) {
	for (size_t i = 0; i < num_points; i++) {
		float res = atan2f(ys[i], xs[i]);

		//if (ys[i] < 0.0f)
			//out[i] += 2 * M_PI;
		out[i] = (ys[i] < 0.0f ? 2 * M_PI : 0) + res;

	if(i == num_points - 1)
		printf("Baseline atan2(%.2f/%.2f)=%.2f\n", ys[i], xs[i], out[i]*180/M_PI);
	}
}

inline float atan_scalar_approximation(float x) {
	float a1 = 0.99997726f;
	float a3 = -0.33262347f;
	float a5 = 0.19354346f;
	float a7 = -0.11643287f;
	float a9 = 0.05265332f;
	float a11 = -0.01172120f;

	float x_sq = x*x;
	return
		x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

inline float atan_scalar_approximation_global(float x) {
	float x_sq = x*x;
	return
		x * (A1 + x_sq * (A3 + x_sq * (A5 + x_sq * (A7 + x_sq * (A9 + x_sq * A11)))));
}

void atan2_auto_1(size_t num_points, const float* ys, const float* xs, float* out) {
	float pi = M_PI;
	float pi_2 = M_PI_2;

	for (size_t i = 0; i < num_points; i++) {
		// Ensure input is in [-1, +1]
		float y = ys[i];
		float x = xs[i];
		bool swap = fabs(x) < fabs(y);
		float atan_input = (swap ? x : y) / (swap ? y : x);

		// Approximate atan
		float res = atan_scalar_approximation(atan_input);

		// If swapped, adjust atan output
		res = swap ? (atan_input >= 0.0f ? pi_2 : -pi_2) - res : res;

		// Adjust quadrants
		if (x >= 0.0f && y >= 0.0f) {}                     // 1st quadrant
		else if (x <  0.0f && y >= 0.0f) { res = pi + res; } // 2nd quadrant
		else if (x <  0.0f && y <  0.0f) { res = pi + res; } // 3rd quadrant
		else if (x >= 0.0f && y < 0.0f) { res = 2 * pi + res; }                     // 4th quadrant

																// Store result
		out[i] = res;
		//printf("atan2 auto 1 (%.2f/%.2f)=%.2f \n", ys[i], xs[i], out[i] * 180 / M_PI);
	}
}

void atan2_auto_2(size_t num_points, const float* ys, const float* xs, float* out) {
	float pi = M_PI;
	float pi_2 = M_PI_2;

	for (size_t i = 0; i < num_points; i++) {
		// Ensure input is in [-1, +1]
		float y = ys[i];
		float x = xs[i];
		bool swap = fabs(x) < fabs(y);
		float atan_input = (swap ? x : y) / (swap ? y : x);

		// Approximate atan
		float res = atan_scalar_approximation(atan_input);
		// If swapped, adjust atan output
		res = swap ? (atan_input >= 0.0f ? pi_2 : -pi_2) - res : res;
		// Adjust quadrants
		if (x >= 0.0f) { res = (y >= 0.0f ? 0 : 2 * pi) + res; }                     // 1st quadrant
		//else if (x <  0.0f) { res = pi + res; }
		else { res = pi + res; }
		out[i] = res;
		if(i == num_points -1)
			printf("atan2 auto 2(%.2f/%.2f)=%.2f \n", ys[i], xs[i], out[i] * 180 / M_PI);
	}
}

void atan2_auto_2_onebyone(const float ys, const float xs, float out) {
	float pi = M_PI;
	float pi_2 = M_PI_2;

	//for (size_t i = 0; i < num_points; i++) {
		// Ensure input is in [-1, +1]
		float y = ys;
		float x = xs;
		bool swap = fabs(x) < fabs(y);
		float atan_input = (swap ? x : y) / (swap ? y : x);

		// Approximate atan
		float res = atan_scalar_approximation_global(atan_input);
		// If swapped, adjust atan output
		res = swap ? (atan_input >= 0.0f ? pi_2 : -pi_2) - res : res;
		// Adjust quadrants
		if (x >= 0.0f) { res = (y >= 0.0f ? 0 : 2 * pi) + res; }                     // 1st quadrant
																					 //else if (x <  0.0f) { res = pi + res; }
		else { res = pi + res; }
		out = res;
		//printf("atan2 auto 2 one by one (%.2f/%.2f)=%.2f \n", y, x, out * 180 / M_PI);
	//}
}

double GetATANDeg(double dY, double dX)
{
	if (dX == 0 && dY>0)
		return 90.0;
	if (dX == 0 && dY<0)
		return 270.0;
	if (dY == 0 && dX>0)
		return 0.0;
	if (dY == 0 && dX<0)
		return 180.0;
	double dVal = abs((dY) / (dX));
	double dPhaseDeg = 999.0;
	int iIterStart = 800;
	for (int j = 100; j<899; j += 100)
	{
		if (dVal<ATAN_LOOKUP[j])
		{
			iIterStart = j - 100;
			break;
		}
	}
	//for(int i=iIterStart;i<iIterStart+100;i+=10) // 081013
	bool bFound = false;
	for (int i = iIterStart; i<899; i += 10)
	{
		if (dVal<ATAN_LOOKUP[i])
		{
			iIterStart = i - 10;
			bFound = true;
			break;
		}
	}
	if (!bFound && (iIterStart == 800)) // Special case to cater last 10 position in the look up table
		iIterStart = 890;
	for (int i = iIterStart; i<899; i++)
	{
		if (ATAN_LOOKUP[i]>dVal)
		{
			if (ATAN_LOOKUP[i] - dVal <= dVal - ATAN_LOOKUP[i - 1])
			{
				dPhaseDeg = ((double)i + 1.0)*0.1;
				break;
			}
			else
			{
				dPhaseDeg = (double)i*0.1;
				break;
			}
		}
	}
	if (dPhaseDeg == 999.0)
		dPhaseDeg = 89.95;
	if (dY>0 && dX>0) // no changes -- first quad
		return dPhaseDeg;
	else if (dY>0 && dX<0) // 2nd quad
		return 180.0 - dPhaseDeg;
	else if (dY<0 && dX<0) // 3rd quad
		return dPhaseDeg + 180.0;
	else // 4th quad
		return 360.0 - dPhaseDeg;
}
