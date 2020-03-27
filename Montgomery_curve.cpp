#include <cstring>
#include "Montgomery_curve.hpp"

Point::Point() {
    mpz_inits(X, Z, NULL);
}

Point::~Point() {
    mpz_clears(X, Z, NULL);
}

void Point::set_as_P() {
    //Так как предполагаем, что Z = 1, то X есть x - координата точки P в афинной форме.
    //Переведенная из координат для кривой Эдвардса в координаты для кривой Монтгомери x - координата точки P:
    mpz_set_str(X, "CBB8F5EBD80486B923EBFB17E5464173144CAC7B0447717B0EA8DE20545A6A23", 16);
    //Соответственно, Z = 1:
    mpz_set_ui(Z, 1ul);
}

void Point::set_as_zero() {
    //Нулевая точка - это (1:0), её удвоение дает снова (1:0),
    //переход от n (n=0 в данном случае) к 2n+1 (2*0+1=1) для любой
    //ненулевой точки Q дает Q, однако к (1:0) неприменима операция перехода от
    //n (n=0 в данном случае) к 2n+1, если следующая точка (2*0+1=1) снова полагается равной (1:0)
    mpz_set_ui(X, 1ul);
    mpz_set_ui(Z, 0ul);
}

Point& Point::operator=(const Point &other) {
    mpz_set(X, other.X);
    mpz_set(Z, other.Z);
    return *this;
}

MontgomeryCurve::MontgomeryCurve() {
    mpz_init_set_str(p, "00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 16);
    //Переведенные из формы Эдвардса в форму Монтгомери коэффициенты:
    mpz_init_set_str(A, "ABCD6AB42CF78BD83F256AE2E7089E30F31637086D4E41FD4ECF952F2B6C6E86", 16);
    mpz_init_set_str(B, "ABCD6AB42CF78BD83F256AE2E7089E30F31637086D4E41FD4ECF952F2B6C6E88", 16);
    mpz_init_set_str(m, "01000000000000000000000000000000003F63377F21ED98D70456BD55B0D8319C", 16);
    mpz_init_set_str(q, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 16);
    //Переведенные из координат для кривой Эдвардса в координаты для кривой Монтгомери координаты точки P:
    mpz_init_set_str(x, "CBB8F5EBD80486B923EBFB17E5464173144CAC7B0447717B0EA8DE20545A6A23", 16);
    mpz_init_set_str(y, "370E3A4D38005921EF122701D68F401C8B685C09767BA6448AF94C29DF1AA555", 16);
    //Посчитанное заранее значение (A-2)/4 mod p:
    mpz_init_set_str(C, "2AF35AAD0B3DE2F60FC95AB8B9C2278C3CC58DC21B53907F53B3E54BCADB1BA1", 16);
    //Проинициализируем остальные объекты:
    P.set_as_P();
}

MontgomeryCurve::~MontgomeryCurve() {
    mpz_clears(p, A, B, m, q, x, y, C, NULL);
}

void MontgomeryCurve::get_q(mpz_t rop) const {
    mpz_set(rop, q);
}

void MontgomeryCurve::doubling(Point &nP) const {
    //Формулы для перехода от [n]P = (Xn:Zn) к [2n]P = (X2n:Z2n):
    //C = (A-2)/4 mod p
    //sqrXpZ = (Xn + Zn)^2
    //sqrXmZ = (Xn - Zn)^2
    //sqrDif = sqrXpZ - sqrXmZ
    //X2n = sqrXmZ * sqrXpZ mod p
    //Z2n = sqrDif * (sqrXpZ + C * sqrDif) mod p

    mpz_t sqrXpZ, sqrXmZ, sqrDif; //Временные объекты
    mpz_inits(sqrXpZ, sqrXmZ, sqrDif, NULL);

    mpz_add(sqrXpZ, nP.X, nP.Z); //sqrXpZ = Xn + Zn
    mpz_pow_ui(sqrXpZ, sqrXpZ, 2ul); //sqrXpZ = sqrXpZ^2 = (Xn + Zn)^2

    mpz_sub(sqrXmZ, nP.X, nP.Z); //sqrXmZ = Xn - Zn
    mpz_pow_ui(sqrXmZ, sqrXmZ, 2ul); //sqrXmZ = sqrXmZ^2 = (Xn - Zn)^2

    mpz_sub(sqrDif, sqrXpZ, sqrXmZ); //sqrDif = sqrXpZ - sqrXmZ

    mpz_mul(nP.X, sqrXmZ, sqrXpZ); //X = sqrXmZ * sqrXpZ
    mpz_mod(nP.X, nP.X, p); //X = X2n = sqrXmZ * sqrXpZ mod p

    mpz_addmul(sqrXpZ, C, sqrDif); //sqrXpZ = sqrXpZ + C * sqrDif
    mpz_mul(nP.Z, sqrDif, sqrXpZ); //Z = sqrDif * (sqrXpZ + C * sqrDif)
    mpz_mod(nP.Z, nP.Z, p); //Z = Z2n = sqrDif * (sqrXpZ + C * sqrDif) mod p

    mpz_clears(sqrXpZ, sqrXmZ, sqrDif, NULL);
}

void MontgomeryCurve::doubling_and_inc(Point &nP, const Point &nPnext, const Point &firstP) const {
    //Формулы для перехода от [n]P и [n+1]P к [2n+1]P:
    //subAddMul = (Xn - Zn)*( X(n+1) + Z(n+1) )
    //addSubMul = (Xn + Zn)*( X(n+1) - Z(n+1) )
    //X(2n+1) = (subAddMul + addSubMul)^2 * Z1 mod p
    //Z(2n+1) = (subAddMul - addSubMul)^2 * X1 mod p

    mpz_t subAddMul, addSubMul, tmp;
    mpz_inits(subAddMul, addSubMul, tmp, NULL);

    mpz_sub(subAddMul, nP.X, nP.Z); //subAddMul = Xn - Zn
    mpz_add(tmp, nPnext.X, nPnext.Z); //tmp = X(n+1) + Z(n+1)
    mpz_mul(subAddMul, subAddMul, tmp); //subAddMul = subAddMul*tmp = (Xn - Zn)*( X(n+1) + Z(n+1) )

    mpz_add(addSubMul, nP.X, nP.Z); //addSubMul = Xn + Zn
    mpz_sub(tmp, nPnext.X, nPnext.Z); //tmp = X(n+1) - Z(n+1)
    mpz_mul(addSubMul, addSubMul, tmp); //addSubMul = addSubMul*tmp = (Xn + Zn)*( X(n+1) - Z(n+1) )

    mpz_add(nP.X, subAddMul, addSubMul); //X = subAddMul + addSubMul
    mpz_pow_ui(nP.X, nP.X, 2ul); //X = X^2 = (subAddMul + addSubMul)^2
    mpz_mul(nP.X, nP.X, firstP.Z); //X = X * Z1 = (subAddMul + addSubMul)^2 * Z1
    mpz_mod(nP.X, nP.X, p); //X = X(2n+1) = (subAddMul + addSubMul)^2 * Z1 mod p

    mpz_sub(nP.Z, subAddMul, addSubMul); //Z = subAddMul - addSubMul
    mpz_pow_ui(nP.Z, nP.Z, 2ul); //Z = Z^2 = (subAddMul - addSubMul)^2
    mpz_mul(nP.Z, nP.Z, firstP.X); //Z = Z * X1 = (subAddMul - addSubMul)^2 * X1
    mpz_mod(nP.Z, nP.Z, p);//Z = Z(2n+1) = (subAddMul - addSubMul)^2 * X1 mod p

    mpz_clears(subAddMul, addSubMul, tmp, NULL);
}

int MontgomeryCurve::power(Point &nP, const mpz_t pow) const {
    //Если степень или координаты точки отрицательны
    if (mpz_cmp_ui(pow, 0ul) < 0 || mpz_cmp_ui(nP.X, 0ul) < 0 || mpz_cmp_ui(nP.Z, 0ul) < 0)
        return -1;
    //Если точка не принадлежит кривой:
    if(point_belonging(nP) == -1) {
        return -1;
    }
    //Если степень равна нулю или задана нулевая точка (1:0)~(r:0)
    if (mpz_cmp_ui(pow, 0ul) == 0 || (mpz_cmp_ui(nP.X, 1ul) >= 0 && mpz_cmp_ui(nP.Z, 0ul) == 0)) {
        nP.set_as_zero();
        return 0;
    }

    char *strPow = mpz_get_str(NULL, 2, pow); //Двоичное разложение степени
    size_t powLength = strlen(strPow); //Число двоичных разрядов степени
    Point firstP, nPnext, tmp;
    firstP = nP;
    nPnext = nP;
    doubling(nPnext); //nPnext = [2]P

    //Пропускаем старший бит степени, так как уже имеем точку первой степени P:
    for (size_t i = 1; i < powLength; ++i)
        if (strPow[i] == '0') {
            tmp = nP;
            doubling_and_inc(tmp, nPnext, firstP);
            nPnext = tmp;
            doubling(nP);
        } else {
            doubling_and_inc(nP, nPnext, firstP);
            doubling(nPnext);
        }

    mpz_t inverse;
    mpz_init(inverse);
    //Если Zpow != 0, то можем найти обратный элемент (так как модуль p - простой) и перейти к x = Xpow*Zpow^(-1):
    if(mpz_cmp_ui(nP.Z, 0ul) != 0) {
        //И все же проверим результат работы mpz_invert:
        if (!mpz_invert(inverse, nP.Z, p)) {
            mpz_clear(inverse);
            return 0;
        }
        mpz_mul(nP.X, nP.X, inverse);
        mpz_mod(nP.X, nP.X, p);
        //По идее, можно сразу присвоить Zpow = 1, но сделаем все "по-честному" только лишь с целью проверки
        //правильности вычислений:
        mpz_mul(nP.Z, nP.Z, inverse);
        mpz_mod(nP.Z, nP.Z, p);
    } else
        if (mpz_cmp_ui(nP.X, 0ul) != 0) {
            //Если Zpow = 0, и Xpow != 0, то можем перейти к Xpow = 1, так как существует обратный к Xpow:
            //Опять же, сделаем все "по-честному":
            if (!mpz_invert(inverse, nP.X, p)) {
                mpz_clear(inverse);
                return 0;
            }
            mpz_mul(nP.X, nP.X, inverse);
            mpz_mod(nP.X, nP.X, p);
        }
    mpz_clear(inverse);
    return 0;
}

int MontgomeryCurve::point_belonging(const Point &point) const {
    mpz_t a, tmp;
    //a есть правая часть уравнения кривой, т.е. a = x^3+A*x^2+x
    mpz_inits(a, tmp, NULL);
    mpz_pow_ui(tmp, point.X, 2ul); //tmp = x^2
    mpz_set(a, point.X); //a = x
    mpz_addmul(a, A, tmp); //a = a + A*x^2 = x + A*x^2
    mpz_addmul(a, point.X, tmp); //a = a + x*x^2 = x + A*x^2 + x^3
    mpz_invert(tmp, B, p); //tmp = B^(-1) mod p
    mpz_mul(a, a, tmp); //a = a*B^(-1)
    mpz_mod(a, a, p); //a = a mod p
    int res = mpz_legendre(a, p); //res = -1 или 1, так как a всегда меньше p.
    //Нельзя сделать было сразу return mpz_legendre(), потому что нужно отчистить память:
    mpz_clears(a, tmp, NULL);
    return res;
}