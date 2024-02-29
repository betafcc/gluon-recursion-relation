watch() {
    local file="$1"
    local std="${2:-c++17}"

    ls ./cc/*.cc |
        entr -r -s "clear && g++ -std=${std} ./cc/${file}.cc -o ./build/${file} && echo '\n======\n' && ./build/${file}"
}
