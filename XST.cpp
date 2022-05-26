#include "XST.h"


XST::XST()
{
}

XST::XST(const XST& oldXST)
{
	this->step = oldXST.step;
	this->a[0] = oldXST.a[0];
	this->a[1] = oldXST.a[1];
	this->a[2] = oldXST.a[2];
	this->b[0] = oldXST.b[0];
	this->b[1] = oldXST.b[1];
	this->b[2] = oldXST.b[2];
	this->c[0] = oldXST.c[0];
	this->c[1] = oldXST.c[1];
	this->c[2] = oldXST.c[2];
	this->o[0] = oldXST.o[0];
	this->o[1] = oldXST.o[1];
	this->o[2] = oldXST.o[2];
}


XST::~XST()
{
}

void XST::getStep(int& Step)
{
	Step = this->step;
}

void XST::getA(float* aOut)
{
	aOut[0] = this->a[0];
	aOut[1] = this->a[1];
	aOut[2] = this->a[2];
}

void XST::getB(float* bOut)
{
	bOut[0] = this->b[0];
	bOut[1] = this->b[1];
	bOut[2] = this->b[2];
}

void XST::getC(float* cOut)
{
	cOut[0] = this->c[0];
	cOut[1] = this->c[1];
	cOut[2] = this->c[2];
}

void XST::getO(float* oOut)
{
	oOut[0] = this->o[0];
	oOut[1] = this->o[1];
	oOut[2] = this->o[2];
}

void XST::setStep(const int& Step)
{
	this->step = Step;
}

void XST::setA(const float* aIn)
{
	this->a[0] = aIn[0];
	this->a[1] = aIn[1];
	this->a[2] = aIn[2];
}

void XST::setB(const float* bIn)
{
	this->b[0] = bIn[0];
	this->b[1] = bIn[1];
	this->b[2] = bIn[2];
}

void XST::setC(const float* cIn)
{
	this->c[0] = cIn[0];
	this->c[1] = cIn[1];
	this->c[2] = cIn[2];
}

void XST::setO(const float* oIn)
{
	this->o[0] = oIn[0];
	this->o[1] = oIn[1];
	this->o[2] = oIn[2];
}