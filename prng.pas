{ generatory psevdosluchaynyh chisel: neyman, modificirovannyy neyman, lemer, kongruentnyy }
{ tekst tolko ascii, bez SysUtils i bez dinamicheskih massivov: sobiraetsya v lyubom fpc }
program Prng;

const
  NeumannDigits = 4;
  MaxVals = 16;
  MaxMod = 4096;

type
  TStrArr = array[0..MaxVals - 1] of string;

{ oshibka s ponyatnym tekstom i vyhod }
procedure Fail(msg: string);
begin
  WriteLn('error: ', msg);
  Halt(1);
end;

{ celoe v stroku }
function IntStr(v: Int64): string;
var
  s: string;
begin
  Str(v:0, s);
  IntStr := s;
end;

{ celoe, dopolnennoe nulyami sleva do w znakov }
function Pad(v: Int64; w: Integer): string;
var
  s: string;
begin
  s := IntStr(v);
  while Length(s) < w do
    s := '0' + s;
  Pad := s;
end;

{ celoe v stroku vida 0,xxxx }
function MakeFrac(v: Int64; d: Integer): string;
begin
  MakeFrac := '0,' + Pad(v, d);
end;

{ d srednih cifr proizvedeniya, dopolnennogo nulyami do 2*d znakov }
function Middle(product: Int64; d: Integer): Int64;
var
  s: string;
  r: Int64;
  code: Integer;
begin
  s := Pad(product, 2 * d);
  if Length(s) <> 2 * d then
    Fail('product ' + IntStr(product) + ' does not fit in ' + IntStr(2 * d) + ' digits');
  { rezhem seredinu: otbrasyvaem d div 2 cifr sleva i sprava }
  Val(Copy(s, d div 2 + 1, d), r, code);
  Middle := r;
end;

{ metod serednyh kvadratov: R(i+1) = seredina(R(i)^2) }
procedure Neumann(r0: Int64; n: Integer; var vals: TStrArr);
var
  i: Integer;
  cur, sq: Int64;
begin
  if n + 1 > MaxVals then
    Fail('too many steps');
  WriteLn('zadacha 1: metod Neymana');
  WriteLn('  R0 = ', MakeFrac(r0, NeumannDigits));
  vals[0] := MakeFrac(r0, NeumannDigits);
  cur := r0;
  for i := 1 to n do
  begin
    { vozvodim v kvadrat i beryom srednie cifry }
    sq := cur * cur;
    cur := Middle(sq, NeumannDigits);
    WriteLn('  R', i - 1, '^2 = ', Pad(sq, 2 * NeumannDigits), ' -> R', i, ' = ',
      MakeFrac(cur, NeumannDigits));
    vals[i] := MakeFrac(cur, NeumannDigits);
  end;
end;

{ modificirovannyy neyman: R(i+1) = seredina(R(i-1) * R(i)) }
procedure NeumannModified(r0, r1: Int64; n: Integer; var vals: TStrArr);
var
  i: Integer;
  prev, cur, prod, nxt: Int64;
begin
  if n + 2 > MaxVals then
    Fail('too many steps');
  WriteLn('zadacha 2: modificirovannyy neyman');
  WriteLn('  R0 = ', MakeFrac(r0, NeumannDigits));
  WriteLn('  R1 = ', MakeFrac(r1, NeumannDigits));
  vals[0] := MakeFrac(r0, NeumannDigits);
  vals[1] := MakeFrac(r1, NeumannDigits);
  prev := r0;
  cur := r1;
  for i := 2 to n + 1 do
  begin
    { peremnozhaem dva predyduschih chisla i beryom seredinu }
    prod := prev * cur;
    nxt := Middle(prod, NeumannDigits);
    WriteLn('  R', i - 2, '*R', i - 1, ' = ', Pad(prod, 2 * NeumannDigits), ' -> R', i,
      ' = ', MakeFrac(nxt, NeumannDigits));
    vals[i] := MakeFrac(nxt, NeumannDigits);
    prev := cur;
    cur := nxt;
  end;
end;

{ algoritm Lemera y = frac(g*x); vsyo v celyh, bez real }
procedure Lehmer(r0, g: Int64; n: Integer; var vals: TStrArr);
const
  Digits = 3;
  Scale = 1000;
var
  i: Integer;
  cur, full: Int64;
begin
  if n + 1 > MaxVals then
    Fail('too many steps');
  WriteLn('zadacha 3: algoritm Lemera {g*x}');
  WriteLn('  R0 = ', MakeFrac(r0, Digits));
  vals[0] := MakeFrac(r0, Digits);
  cur := r0;
  for i := 1 to n do
  begin
    { g*R v celyh: celaya chast uhodit, ostayotsya drobnaya }
    full := g * cur;
    cur := full mod Scale;
    WriteLn('  g*R', i - 1, ' = ', full div Scale, ',', Pad(cur, Digits), ' -> R', i,
      ' = ', MakeFrac(cur, Digits));
    vals[i] := MakeFrac(cur, Digits);
  end;
end;

{ x/m s okrugleniem half-up (tak schitayut na bumage), tri znaka }
function Ratio(x, m: Int64): string;
var
  r: Int64;
begin
  r := (x * 2000 + m) div (2 * m);
  Ratio := IntStr(r div 1000) + ',' + Pad(r mod 1000, 3);
end;

{ multiplikativnyy kongruentnyy metod X(i+1) = a*X(i) mod m;
  period = -1, esli cikl ne zamknulsya }
procedure Congruential(a, m, x0: Int64; n: Integer; var vals: TStrArr; var period: Integer);
var
  i, steps: Integer;
  seen: array[0..MaxMod - 1] of Integer;
  cur, prev, probe: Int64;
begin
  if (m <= 0) or (m > MaxMod) then
    Fail('modulus m must be in 1..' + IntStr(MaxMod) + ', got ' + IntStr(m));
  if (x0 < 0) or (x0 >= m) then
    Fail('X0 must be in 0..' + IntStr(m - 1) + ', got ' + IntStr(x0));
  if n > MaxVals then
    Fail('too many steps');

  { seen[x] = nomer shaga, na kotorom x uzhe vstrechalsya, -1 esli net }
  for i := 0 to m - 1 do
    seen[i] := -1;
  seen[x0] := 0;
  period := -1;

  WriteLn('zadacha 4: kongruentnyy metod a=', a, ', m=', m);
  WriteLn('  X0 = ', x0);
  cur := x0;
  for i := 1 to n do
  begin
    { delim a*X na m: chastnoe i ostatok (ostatok i est sleduyuschiy X) }
    prev := cur;
    cur := (a * prev) mod m;
    vals[i - 1] := Ratio(cur, m);
    WriteLn('  a*X', i - 1, '/m = ', a, '*', prev, '/', m, ' = ', (a * prev) div m, '*', m,
      '+', cur, ' -> X', i, ' = ', cur, ', R', i, ' = ', vals[i - 1]);

    { fiksiruem vozvrat v uzhe videnoe sostoyanie }
    if (period < 0) and (seen[cur] >= 0) then
      period := i - seen[cur];
    if seen[cur] < 0 then
      seen[cur] := i;
  end;

  { esli za n shagov cikl ne zamknulsya, dogonyaem do zamykaniya (period ne bolshe m) }
  probe := cur;
  steps := n;
  while (period < 0) and (steps < m + 1) do
  begin
    probe := (a * probe) mod m;
    Inc(steps);
    if seen[probe] >= 0 then
      period := steps - seen[probe]
    else
      seen[probe] := steps;
  end;
end;

{ sverka s etalonom iz metodichki: znacheniya cherez probel odnoy strokoy }
procedure Check(lbl: string; var got: TStrArr; count: Integer; expected: string);
var
  i: Integer;
  joined: string;
begin
  joined := '';
  for i := 0 to count - 1 do
  begin
    if i > 0 then
      joined := joined + ' ';
    joined := joined + got[i];
  end;
  if joined <> expected then
    Fail(lbl + ': expected [' + expected + '], got [' + joined + ']');
end;

var
  t1, t2, t3, t4: TStrArr;
  period: Integer;

begin
  { zadacha 1: neyman, N=6, R0=0,5836 }
  Neumann(5836, 6, t1);
  Check('zadacha 1', t1, 7, '0,5836 0,0588 0,3457 0,9508 0,4020 0,1604 0,5728');
  WriteLn;

  { zadacha 2: modificirovannyy neyman, N=6, R0=0,5836, R1=0,2176 }
  NeumannModified(5836, 2176, 6, t2);
  Check('zadacha 2', t2, 8, '0,5836 0,2176 0,6991 0,2124 0,8488 0,0285 0,4190 0,1941');
  WriteLn;

  { zadacha 3: lemer frac(g*x), N=5, R0=0,585, g=927 }
  Lehmer(585, 927, 5, t3);
  Check('zadacha 3', t3, 6, '0,585 0,295 0,465 0,055 0,985 0,095');
  WriteLn;

  { zadacha 4: kongruentnyy metod, a=265, m=129, X0=122, N=5 }
  Congruential(265, 129, 122, 5, t4, period);
  Check('zadacha 4', t4, 5, '0,620 0,341 0,388 0,713 0,992');
  WriteLn;

  { period kongruentnogo metoda vyvodim otdelno }
  WriteLn('period zadachi 4: ', period);
end.
