#ifndef XST_H
#define XST_H
#include <fstream>
using namespace std;
class XST
{
public:
	XST();
	XST(const XST& oldXST);
	~XST();
	void getStep(int& Step);
	void getA(float* aOut);
	void getB(float* bOut);
	void getC(float* cOut);
	void getO(float* oOut);
	void setStep(const int& Step);
	void setA(const float* aIn);
	void setB(const float* bIn);
	void setC(const float* cIn);
	void setO(const float* oIn);

protected:
	int step;
	float a[3], b[3], c[3], o[3];
};

#endif

