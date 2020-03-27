#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "Montgomery_curve.hpp"

void point_input(Point &point, const MontgomeryCurve &curve) {
    char *input = NULL; //Указатель на введенную строку
    size_t cap; //Выделенная емкость под input
    ssize_t size; //Количество записанных в input символов
    printf("Выбор точки:\n");
    begin:
    printf("Если вы желаете использовать точку P (порождающий элемент подгруппы порядка q), введите 'p'\n");
    printf("Если вы желаете задать точку самостоятельно, введите x-координату точки (неотрицательное число)\n");
    printf("Возможные префиксы систем счисления: 0X (0x) - 16-ичная, 0B (0b) - 2-ичная, 0 - 8-ичная, иначе - 10-ичная\n");

    size = getline(&input, &cap, stdin);
    //Удалим \n из конца ввода:
    input[size-1] = '\0';
    //Если введен только символ p + \n:
    if (size == 2 && input[0] == 'p') {
        point.set_as_P();
    } else {
        //Иначе проверим, введено ли корретное неотрицательное число:
        //(set_str возвращает -1 в случае некорректности строки, 0 - в случае успеха):
        if (input[0] == '-' || mpz_set_str(point.X, input, 0)) {
            printf("Ошибка ввода\n");
            goto begin;
        }
        //X введен, зададим Z = 1:
        mpz_set_ui(point.Z, 1ul);
        if(curve.point_belonging(point) == -1) {
            printf("Точка не принадлежит кривой\n");
            goto begin;
        }
    }
    free(input);
}

void power_input(mpz_t pow, mpz_t q, mpz_t p, char *powLabel) {
    char *input = NULL; //Указатель на введенную строку
    size_t cap; //Выделенная емкость под input
    ssize_t size; //Количество записанных в input символов
    printf("Выбор степени:\n");
    begin:
    printf("Если вы хотите возвести точку в степень q, введите 'q', иначе:\n");
    printf("Введите 'q+', обозначающий сложение числа с q, 'q-' - вычитание из q, 'q*' - умножение q на число\n");
    printf("После выбранной команды введите необходимое положительное число, например: 'q+2' - возведение точки в степень (q+2)\n");
    printf("Если вы желаете задать степень самостоятельно, просто введите нужное неотрицательное число\n");
    printf("Возможные префиксы систем счисления: 0X (0x) - 16-ичная, 0B (0b) - 2-ичная, 0 - 8-ичная, иначе - 10-ичная\n");

    size = getline(&input, &cap, stdin);
    //Удалим \n из конца ввода:
    input[size-1] = '\0';
    //Если не введен символ q:
    if(input[0] != 'q') {
        //Тогда проверим, введено ли корретное неотрицательное число:
        //(set_str возвращает -1 в случае некорректности строки, 0 - в случае успеха):
        if (input[0] == '-' || mpz_set_str(p, input, 0)) {
            printf("Ошибка ввода\n");
            goto begin;
        }
        //В данном случае p уже и есть необходимая степень, скопируем значение в pow,
        //так как финальный вывод имеет один шаблон для всех случаев:
        mpz_set(pow, p);
        strcpy(powLabel, "(p)");
    } else {
        //Если первым введен символ q и \n:
        if (size == 2) {
            //По сути, p в данном случае не нужно, но зададим его, так как финальный вывод имеет один шаблон для всех случаев:
            mpz_set(pow, q);
            mpz_set_ui(p, 0);
            strcpy(powLabel, "(q)");
            free(input);
            return;
        }
        //Иначе проверим, введены ли символы +,-,* и корретное неотрицательное число:
        if((input[1] != '+' && input[1] != '-' && input[1] != '*') || input[2] == '-' || mpz_set_str(p, input+2, 0)) {
            printf("Ошибка ввода\n");
            goto begin;
        }
        switch(input[1]) {
            case '+':
                mpz_add(pow, q, p);
                strcpy(powLabel, "(q+p)");
                break;
            case '-':
                mpz_sub(pow, q, p);
                strcpy(powLabel, "(q-p)");
                break;
            case '*':
                mpz_mul(pow, q, p);
                strcpy(powLabel, "(q*p)");
        }
    }
    free(input);
}

int main() {
    MontgomeryCurve curve;
    Point point;
    mpz_t pow, q, p; //pow - степень для возведения в неё точки; q - порядок P; p - число для модификации q
    mpz_inits(pow, q, p, NULL);
    curve.get_q(q);
    char powLabel[6] = {}; //строка для наглядного отображения при выводе вычисляемой степени
    printf("Данная программа возводит в указанную степень указанную точку на эллиптической кривой в форме Монтгомери");
    printf("\nс заданными коэффициентами A, B и по заданному модулю p.\n\n");
    point_input(point, curve);
    power_input(pow, q, p, powLabel);

    printf("X(1) = ");
    mpz_out_str(stdout, 10, point.X);
    printf("\nZ(1) = ");
    mpz_out_str(stdout, 10, point.Z);
    printf("\n\n");

    if(curve.power(point, pow)) {
        printf("Ошибка вычисления\n");
        return -1;
    }
    printf("p = ");
    mpz_out_str(stdout, 10, p);
    printf("\nX%s = ", powLabel);
    mpz_out_str(stdout, 10, point.X);
    printf("\nZ%s = ", powLabel);
    mpz_out_str(stdout, 10, point.Z);

    mpz_clears(pow, q, p, NULL);
    return 0;
}