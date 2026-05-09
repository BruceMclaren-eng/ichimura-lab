#include <iostream>

int main() {
    int a = 42;

    int* p = &a;
    int& r = a;

    std::cout << "aの値:     " << a <<"\n";
    std::cout << "aのアドレス " << &a << "\n";
    std::cout << "pの値      " << p << "\n";
    std::cout << "pが指す値   "<< *p << "\n";
    std::cout << "rの値      "<< r << "\n";
    
    *p = 100;
    std::cout << "変更後の a: "<< a << "\n";

    r = 200;
    std::cout << "変更後の a: "<< a << "\n";

    return 0;
}