"""генераторы псевдослучайных чисел: нейман, модифицированный нейман, лемер, конгруэнтный.

совместимо с python 3.4+: без f-строк, dataclasses и аннотаций переменных.
"""

from collections import namedtuple
from decimal import ROUND_HALF_UP, Decimal

# ширина дробной части для методов неймана (4 знака)
NEUMANN_DIGITS = 4

# результат одной задачи: заголовок, строки решения и значения R
Task = namedtuple("Task", ["title", "rows", "values"])


def coerce_fraction(value: int, digits: int) -> str:
    """превращает целое в строку вида 0,xxxx с нужным числом знаков."""
    return "0,{0:0{1}d}".format(value, digits)


def make_middle(product: int, digits: int) -> int:
    """берёт `digits` средних цифр у произведения, дополненного нулями до 2*digits знаков."""
    # дополняем нулями слева, чтобы середина считалась одинаково для всех чисел
    padded = "{0:0{1}d}".format(product, 2 * digits)
    if len(padded) != 2 * digits:
        raise ValueError("произведение {0} не влезает в {1} знаков".format(product, 2 * digits))

    # режем середину: отбрасываем digits//2 цифр слева и справа
    cut = digits // 2
    return int(padded[cut:cut + digits])


def resolve_neumann(r0: int, n: int, digits: int = NEUMANN_DIGITS) -> Task:
    """метод серединных квадратов: R(i+1) = середина(R(i)^2)."""
    rows = ["R0 = {0}".format(coerce_fraction(r0, digits))]
    values = [coerce_fraction(r0, digits)]
    current = r0

    for i in range(1, n + 1):
        # возводим в квадрат и выбираем средние цифры
        square = current ** 2
        current = make_middle(square, digits)
        rows.append("R{0}^2 = {1:0{2}d} -> R{3} = {4}".format(
            i - 1, square, 2 * digits, i, coerce_fraction(current, digits)))
        values.append(coerce_fraction(current, digits))

    return Task("метод неймана", rows, values)


def resolve_neumann_modified(r0: int, r1: int, n: int, digits: int = NEUMANN_DIGITS) -> Task:
    """модифицированный нейман: R(i+1) = середина(R(i-1) * R(i))."""
    rows = [
        "R0 = {0}".format(coerce_fraction(r0, digits)),
        "R1 = {0}".format(coerce_fraction(r1, digits)),
    ]
    values = [coerce_fraction(r0, digits), coerce_fraction(r1, digits)]
    prev, current = r0, r1

    for i in range(2, n + 2):
        # перемножаем два предыдущих числа и берём середину
        product = prev * current
        nxt = make_middle(product, digits)
        rows.append("R{0}*R{1} = {2:0{3}d} -> R{4} = {5}".format(
            i - 2, i - 1, product, 2 * digits, i, coerce_fraction(nxt, digits)))
        values.append(coerce_fraction(nxt, digits))
        prev, current = current, nxt

    return Task("модифицированный нейман", rows, values)


def resolve_lehmer_fraction(r0: int, g: int, n: int, digits: int = 3) -> Task:
    """алгоритм лемера y = {g*x}; работаем в целых, чтобы не ловить ошибки float."""
    scale = 10 ** digits
    rows = ["R0 = {0}".format(coerce_fraction(r0, digits))]
    values = [coerce_fraction(r0, digits)]
    current = r0

    for i in range(1, n + 1):
        # g*R в целых: целая часть уходит, остаётся дробная
        full = g * current
        current = full % scale
        rows.append("g*R{0} = {1},{2:0{3}d} -> R{4} = {5}".format(
            i - 1, full // scale, full % scale, digits, i, coerce_fraction(current, digits)))
        values.append(coerce_fraction(current, digits))

    return Task("алгоритм лемера {g*x}", rows, values)


def make_ratio(x: int, m: int, digits: int = 3) -> str:
    """x/m с округлением half-up (так считают на бумаге), а не банковским."""
    quant = Decimal(1).scaleb(-digits)
    value = (Decimal(x) / Decimal(m)).quantize(quant, rounding=ROUND_HALF_UP)
    return str(value).replace(".", ",")


def resolve_congruential(a: int, m: int, x0: int, n: int) -> tuple:
    """мультипликативный конгруэнтный метод X(i+1) = a*X(i) mod m.

    возвращает пару (Task, период); период None, если цикл не замкнулся.
    """
    if m <= 0:
        raise ValueError("модуль m должен быть положительным, получили {0}".format(m))

    rows = ["X0 = {0}".format(x0)]
    values = []
    seen = {x0: 0}
    period = None
    current = x0

    for i in range(1, n + 1):
        # делим a*X на m: частное и остаток (остаток и есть следующий X)
        prev = current
        quotient, current = divmod(a * prev, m)
        ratio = make_ratio(current, m)
        rows.append("a*X{0}/m = {1}*{2}/{3} = {4}*{3}+{5} -> X{6} = {5}, R{6} = {7}".format(
            i - 1, a, prev, m, quotient, current, i, ratio))
        values.append(ratio)

        # фиксируем возврат в уже виденное состояние
        if period is None and current in seen:
            period = i - seen[current]
        seen.setdefault(current, i)

    # если за n шагов цикл не замкнулся, догоняем до замыкания (период не больше m)
    probe = current
    steps = n
    while period is None and steps < m + 1:
        probe = (a * probe) % m
        steps += 1
        if probe in seen:
            period = steps - seen[probe]
            break
        seen[probe] = steps

    title = "конгруэнтный метод a={0}, m={1}".format(a, m)
    return Task(title, rows, values), period


def check_expected(label: str, got: list, expected: list) -> None:
    """сверка с эталоном из методички, падаем с понятным текстом."""
    if got != expected:
        raise AssertionError("{0}: ожидали {1}, получили {2}".format(label, expected, got))


def main() -> None:
    """прогоняет все четыре задачи и сверяет результат с методичкой."""
    # задача 1: нейман, N=6, R0=0,5836
    t1 = resolve_neumann(5836, 6)
    check_expected(
        "задача 1",
        t1.values,
        ["0,5836", "0,0588", "0,3457", "0,9508", "0,4020", "0,1604", "0,5728"],
    )

    # задача 2: модифицированный нейман, N=6, R0=0,5836, R1=0,2176
    t2 = resolve_neumann_modified(5836, 2176, 6)
    check_expected(
        "задача 2",
        t2.values,
        ["0,5836", "0,2176", "0,6991", "0,2124", "0,8488", "0,0285", "0,4190", "0,1941"],
    )

    # задача 3: лемер {g*x}, N=5, R0=0,585, g=927
    t3 = resolve_lehmer_fraction(585, 927, 5)
    check_expected("задача 3", t3.values, ["0,585", "0,295", "0,465", "0,055", "0,985", "0,095"])

    # задача 4: конгруэнтный метод, a=265, m=129, X0=122, N=5
    t4, period = resolve_congruential(265, 129, 122, 5)
    check_expected("задача 4", t4.values, ["0,620", "0,341", "0,388", "0,713", "0,992"])

    # печатаем решения по порядку
    for idx, task in enumerate((t1, t2, t3, t4), start=1):
        print("задача {0}: {1}".format(idx, task.title))
        for row in task.rows:
            print("  " + row)
        print()

    # период конгруэнтного метода выводим отдельно
    print("период задачи 4: {0}".format(period))


if __name__ == "__main__":
    main()
