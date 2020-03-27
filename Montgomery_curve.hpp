#ifndef MONTGOMERY_ELLIPTIC_CURVE_MONTGOMERY_CURVE_HPP
#define MONTGOMERY_ELLIPTIC_CURVE_MONTGOMERY_CURVE_HPP

#include <gmp.h>

//Структура, описывающая точку (X:Z) в проективном пространстве, где x = X/Z
struct Point {
    mpz_t X;
    mpz_t Z;

    Point();
    ~Point();
    void set_as_P();
    void set_as_zero();
    Point& operator=(const Point &other);
};

//Класс, описывающий эллиптическую кривую в форме Монтгомери и операцию возведения точки в степень на ней:
class MontgomeryCurve {
    mpz_t p; //модуль эллиптической кривой;
    mpz_t A; //коэффициенты уравнения эллиптической кривой
    mpz_t B; //в форме Монтгомери;
    mpz_t m; //порядок группы точек эллиптической кривой;
    mpz_t q; //порядок циклической подгруппы группы точек эллиптической кривой;
    mpz_t x; //координаты точки P (порождающего элемента подгруппы порядка q)
    mpz_t y; //на эллиптической кривой в форме Монтгомери.
    mpz_t C; //Константа (A-2)/4 mod p для удвоения с инкрементом
    Point P; //координаты точки P в проективном пространстве.
    void doubling(Point &point) const;
    void doubling_and_inc(Point &nP, const Point &nPnext, const Point &difP) const;
public:
    MontgomeryCurve();
    ~MontgomeryCurve();
    void get_q(mpz_t rop) const;
    int power(Point &nP, const mpz_t pow) const; //Возвращает 0 в случае успеха, -1 - иначе
    int point_belonging(const Point &point) const;
};

#endif //MONTGOMERY_ELLIPTIC_CURVE_MONTGOMERY_CURVE_HPP