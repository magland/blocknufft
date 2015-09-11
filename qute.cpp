#include "qute.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <chrono>

using namespace std;

//typedef chrono::high_resolution_clock Clock;
class QTimePrivate {
public:
    QTime *q;

    //Clock::time_point m_start_time;

};

QTime::QTime() {
    d=new QTimePrivate;
    d->q=this;
}

QTime::~QTime()
{
    delete d;
}

void QTime::start()
{
    //d->m_start_time=Clock::now();
}

int QTime::elapsed()
{
    //Clock::time_point t2 = Clock::now();
    //long int ms=std::chrono::duration_cast<std::chrono::milliseconds>(t2-d->m_start_time).count();
    //return ms;
}

int qrand()
{
    return rand();
}

class QStringPrivate {
public:
    QString *q;
};

QString::QString() {
    d=new QStringPrivate;
    d->q=this;
}

QString::~QString()
{
    delete d;
}
