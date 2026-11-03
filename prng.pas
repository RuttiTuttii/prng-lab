{ генераторы псевдослучайных чисел: нейман, модифицированный нейман, лемер, конгруэнтный.
  сборка: fpc prng.pas }
program Prng;

{$mode objfpc}{$H+}

uses
  SysUtils;

const
  { ширина дробной части для методов неймана (4 знака) }
  NeumannDigits = 4;

type
  TStrArr = array of string;

{ целое, дополненное нулями слева до w знаков }
function Pad(v: Int64; w: Integer): string;
begin
  Result := IntToStr(v);
  while Length(Result) < w do
    Result := '0' + Result;
end;

{ превращает целое в строку вида 0,xxxx с нужным числом знаков }
function Frac(v: Int64; d: Integer): string;
begin
  Result := '0,' + Pad(v, d);
end;

{ берёт d средних цифр у произведения, дополненного нулями до 2*d знаков }
function Middle(product: Int64; d: Integer): Int64;
var
  s: string;
begin
  s := Pad(product, 2 * d);
  if Length(s) <> 2 * d then
    raise Exception.CreateFmt('произведение %d не влезает в %d знаков', [product, 2 * d]);
  { режем середину: отбрасываем d div 2 цифр слева и справа }
  Result := StrToInt64(Copy(s, d div 2 + 1, d));
end;

{ метод серединных квадратов: R(i+1) = середина(R(i)^2) }
function Neumann(r0: Int64; n: Integer): TStrArr;
var
  i: Integer;
  cur, sq: Int64;
begin
  SetLength(Result, n + 1);
  WriteLn('задача 1: метод неймана');
  WriteLn('  R0 = ', Frac(r0, NeumannDigits));
  Result[0] := Frac(r0, NeumannDigits);
  cur := r0;
  for i := 1 to n do
  begin
    { возводим в квадрат и выбираем средние цифры }
    sq := cur * cur;
    cur := Middle(sq, NeumannDigits);
    WriteLn(Format('  R%d^2 = %s -> R%d = %s',
      [i - 1, Pad(sq, 2 * NeumannDigits), i, Frac(cur, NeumannDigits)]));
    Result[i] := Frac(cur, NeumannDigits);
  end;
end;

{ модифицированный нейман: R(i+1) = середина(R(i-1) * R(i)) }
function NeumannModified(r0, r1: Int64; n: Integer): TStrArr;
var
  i: Integer;
  prev, cur, prod, nxt: Int64;
begin
  SetLength(Result, n + 2);
  WriteLn('задача 2: модифицированный нейман');
  WriteLn('  R0 = ', Frac(r0, NeumannDigits));
  WriteLn('  R1 = ', Frac(r1, NeumannDigits));
  Result[0] := Frac(r0, NeumannDigits);
  Result[1] := Frac(r1, NeumannDigits);
  prev := r0;
  cur := r1;
  for i := 2 to n + 1 do
  begin
    { перемножаем два предыдущих числа и берём середину }
    prod := prev * cur;
    nxt := Middle(prod, NeumannDigits);
    WriteLn(Format('  R%d*R%d = %s -> R%d = %s',
      [i - 2, i - 1, Pad(prod, 2 * NeumannDigits), i, Frac(nxt, NeumannDigits)]));
    Result[i] := Frac(nxt, NeumannDigits);
    prev := cur;
    cur := nxt;
  end;
end;

{ алгоритм лемера y = {g*x}; всё в целых, без float }
function Lehmer(r0, g: Int64; n: Integer): TStrArr;
const
  Digits = 3;
  Scale = 1000;
var
  i: Integer;
  cur, full: Int64;
begin
  SetLength(Result, n + 1);
  WriteLn('задача 3: алгоритм лемера {g*x}');
  WriteLn('  R0 = ', Frac(r0, Digits));
  Result[0] := Frac(r0, Digits);
  cur := r0;
  for i := 1 to n do
  begin
    { g*R в целых: целая часть уходит, остаётся дробная }
    full := g * cur;
    cur := full mod Scale;
    WriteLn(Format('  g*R%d = %d,%s -> R%d = %s',
      [i - 1, full div Scale, Pad(cur, Digits), i, Frac(cur, Digits)]));
    Result[i] := Frac(cur, Digits);
  end;
end;

{ x/m с округлением half-up (так считают на бумаге), три знака }
function Ratio(x, m: Int64): string;
var
  r: Int64;
begin
  r := (x * 2000 + m) div (2 * m);
  Result := IntToStr(r div 1000) + ',' + Pad(r mod 1000, 3);
end;

{ мультипликативный конгруэнтный метод X(i+1) = a*X(i) mod m;
  period = -1, если цикл не замкнулся }
function Congruential(a, m, x0: Int64; n: Integer; out period: Integer): TStrArr;
var
  i, steps: Integer;
  seen: array of Integer;
  cur, prev, probe: Int64;
begin
  if m <= 0 then
    raise Exception.CreateFmt('модуль m должен быть положительным, получили %d', [m]);
  if (x0 < 0) or (x0 >= m) then
    raise Exception.CreateFmt('X0 должен быть в [0, %d), получили %d', [m, x0]);

  { seen[x] = индекс шага, на котором x уже встречался, -1 если нет }
  SetLength(seen, m);
  for i := 0 to m - 1 do
    seen[i] := -1;
  seen[x0] := 0;
  period := -1;

  SetLength(Result, n);
  WriteLn(Format('задача 4: конгруэнтный метод a=%d, m=%d', [a, m]));
  WriteLn('  X0 = ', x0);
  cur := x0;
  for i := 1 to n do
  begin
    { делим a*X на m: частное и остаток (остаток и есть следующий X) }
    prev := cur;
    cur := (a * prev) mod m;
    Result[i - 1] := Ratio(cur, m);
    WriteLn(Format('  a*X%d/m = %d*%d/%d = %d*%d+%d -> X%d = %d, R%d = %s',
      [i - 1, a, prev, m, (a * prev) div m, m, cur, i, cur, i, Result[i - 1]]));

    { фиксируем возврат в уже виденное состояние }
    if (period < 0) and (seen[cur] >= 0) then
      period := i - seen[cur];
    if seen[cur] < 0 then
      seen[cur] := i;
  end;

  { если за n шагов цикл не замкнулся, догоняем до замыкания (период не больше m) }
  probe := cur;
  steps := n;
  while (period < 0) and (steps < m + 1) do
  begin
    probe := (a * probe) mod m;
    Inc(steps);
    if seen[probe] >= 0 then
    begin
      period := steps - seen[probe];
      Break;
    end;
    seen[probe] := steps;
  end;
end;

{ сверка с эталоном из методички, падаем с понятным текстом }
procedure Check(const lbl: string; const got: TStrArr; const expected: array of string);
var
  i: Integer;
begin
  if Length(got) <> Length(expected) then
    raise Exception.CreateFmt('%s: ожидали %d значений, получили %d',
      [lbl, Length(expected), Length(got)]);
  for i := 0 to High(expected) do
    if got[i] <> expected[i] then
      raise Exception.CreateFmt('%s: позиция %d, ожидали %s, получили %s',
        [lbl, i, expected[i], got[i]]);
end;

var
  t1, t2, t3, t4: TStrArr;
  period: Integer;
begin
  { задача 1: нейман, N=6, R0=0,5836 }
  t1 := Neumann(5836, 6);
  Check('задача 1', t1, ['0,5836', '0,0588', '0,3457', '0,9508', '0,4020', '0,1604', '0,5728']);
  WriteLn;

  { задача 2: модифицированный нейман, N=6, R0=0,5836, R1=0,2176 }
  t2 := NeumannModified(5836, 2176, 6);
  Check('задача 2', t2,
    ['0,5836', '0,2176', '0,6991', '0,2124', '0,8488', '0,0285', '0,4190', '0,1941']);
  WriteLn;

  { задача 3: лемер {g*x}, N=5, R0=0,585, g=927 }
  t3 := Lehmer(585, 927, 5);
  Check('задача 3', t3, ['0,585', '0,295', '0,465', '0,055', '0,985', '0,095']);
  WriteLn;

  { задача 4: конгруэнтный метод, a=265, m=129, X0=122, N=5 }
  t4 := Congruential(265, 129, 122, 5, period);
  Check('задача 4', t4, ['0,620', '0,341', '0,388', '0,713', '0,992']);
  WriteLn;

  { период конгруэнтного метода выводим отдельно }
  WriteLn('период задачи 4: ', period);
end.
