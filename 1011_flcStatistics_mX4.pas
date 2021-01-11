//Testunit for maXbox4 scope ! not finished yet
// Study Design for maXbox 4.7.5
//https://github.com/fundamentalslib/fundamentals5/blob/master/Source/Maths/flcMaths.pas
{******************************************************************************}
{                                                                              }
{   Library:          Fundamentals 5.00                                        }
{   File name:        flcStatistics.pas                                        }
{   File version:     5.07                                                     }
{   Description:      Statistic class                                          }
{                                                                              }
{   Copyright:        Copyright (c) 2000-2016, David J Butler                  }
{                     All rights reserved.                                     }
{                     Redistribution and use in source and binary forms, with  }
{                     or without modification, are permitted provided that     }
{                     the following conditions are met:                        }
{                     Redistributions of source code must retain the above     }
{                     copyright notice, this list of conditions and the        }
{                     following disclaimer.                                    }
{                     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND   }
{                     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED          }
{                     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   }
{                     WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A          }
{                     PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL     }
{                     THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    }
{                     INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR             }
{                     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    }
{                     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF     }
{                     USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)         }
{                     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER   }
{                     IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING        }
{                     NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE   }
{                     USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             }
{                     POSSIBILITY OF SUCH DAMAGE.                              }
{                                                                              }
{   Github:           https://github.com/fundamentalslib                       }
{   E-mail:           fundamentals.library at gmail.com                        }
{                                                                              }
{ Revision history:                                                            }
{                                                                              }
{   2000/01/22  0.01  Added TStatistic.                                        }
{   2003/02/16  0.02  Created cStatistics unit from cMaths.                    }
{   2003/03/08  3.03  Revised for Fundamentals 3.                              }
{   2003/03/13  3.04  Skew and Kurtosis. Added documentation.                  }
{   2012/10/26  4.05  Revised for Fundamentals 4.                              }
{   2013/11/17  4.06  Fix bug reported by Steve Eicker.                        }
{   2016/01/17  5.07  Revised for Fundamentals 5.                              }
{                                                                              }
{ Supported compilers:                                                         }
{                                                                              }
{   Delphi XE7 Win32                    5.07  2016/01/17                       }
{   Delphi XE7 Win64                    5.07  2016/01/17                       }
{                                                                              }
{******************************************************************************}

//{$INCLUDE flcMaths.inc}

unit flcStatistics;

interface

(*uses
  { System }
  SysUtils,

  { Fundamentals }
  flcUtils,
  flcMaths;         *)
  

(*  
 {$INCLUDE ..\flcInclude.inc}

{$IFDEF DEBUG}
{$IFDEF TEST}
  {$DEFINE MATHS_TEST}
{$ENDIF}
{$ENDIF}

{$IFDEF ExtendedIsDouble}
  {$DEFINE MFloatIsDouble}
{$ELSE}
  {$DEFINE MFloatIsExtended}
{$ENDIF}   *)


{                                                                              }
{ Exceptions                                                                   }
{                                                                              }
type
  EStatistics = Exception;
  EStatisticsInvalidArgument = {class} EStatistics;
  EStatisticsOverflow = EStatistics;
  
  type MFloat = double;

{                                                                              }
{ Binomial distribution                                                        }
{                                                                              }
{   The binomial distribution gives the probability of obtaining exactly r     }
{   successes in n independent trials, where there are two possible outcomes   }
{   one of which is conventionally called success.                             }
{                                                                              }
{   BinomialCoeff returns the binomial coefficient for the bin(n)-             }
{   distribution.                                                              }
{                                                                              }
function  BinomialCoeff(N, R: Integer): MFloat;



{                                                                              }
{ Normal distribution                                                          }
{                                                                              }
{   The Normal distribution (Gaussian distribution) is a model for values on   }
{   a continuous scale. A normal distribution can be completely described by   }
{   two parameters: mean (m) and variance (s2). It is shown as C ~ N(m, s2).   }
{   The distribution is symmetrical with mean, mode, and median all equal      }
{   at m. In the special case of m = 1 and s2 = 1, it is called the standard   }
{   normal distribution                                                        }
{                                                                              }
{   CumNormal returns the area under the N(u,s) distribution.                  }
{   CumNormal01 returns the area under the N(0,1) distribution.                }
{   InvCummNormal01 returns position on X-axis that gives cummulative area     }
{    of Y0 under the N(0,1) distribution.                                      }
{   InvCummNormal returns position on X-axis that gives cummulative area       }
{     of Y0 under the N(u,s) distribution.                                     }
{                                                                              }
function  erf(x: MFloat): MFloat;
function  erfc(const x: MFloat): MFloat;
function  CummNormal(const u, s, X: MFloat): MFloat;
function  CummNormal01(const X: MFloat): MFloat;
function  InvCummNormal01(Y0: MFloat): MFloat;
function  InvCummNormal(const u, s, Y0: MFloat): MFloat;



{                                                                              }
{ Chi-Squared distribution                                                     }
{                                                                              }
{   The Chi-Squared is a distribution derived from the normal distribution.    }
{   It is the distribution of a sum of squared Normal distributed variables.   }
{   The importance of the Chi-square distribution stems from the fact that it  }
{   describes the distribution of the Variance of a sample taken from a        }
{   Normal distributed population. Chi-squared (C2) is distributed with v      }
{   degrees of freedom with mean = v and variance = 2v.                        }
{                                                                              }
{   CumChiSquare returns the area under the X^2 (Chi-squared) (Chi, Df)        }
{   distribution.                                                              }
{                                                                              }
function  CummChiSquare(const Chi, Df: MFloat): MFloat;



{                                                                              }
{ F distribution                                                               }
{                                                                              }
{   The F-distribution is a continuous probability distribution of the ratio   }
{   of two independent random variables, each having a Chi-squared             }
{   distribution, divided by their respective degrees of freedom. In           }
{   regression analysis, the F-test can be used to test the joint              }
{   significance of all variables of a model.                                  }
{   CumF returns the area under the F (f, Df1, Df2) distribution.              }
{                                                                              }
function  CumF(const f, Df1, Df2: MFloat): MFloat;



{                                                                              }
{ Poison distribution                                                          }
{                                                                              }
{   The Poisson distribution is the probability distribution of the number     }
{   of (rare) occurrences of some random event in an interval of time or       }
{   space. Poisson distribution is used to represent distribution of counts    }
{   like number of defects in a piece of material, customer arrivals,          }
{   insurance claims, incoming telephone calls, or alpha particles emitted.    }
{   A transformation that often changes Poisson data approximately normal is   }
{   the square root.                                                           }
{   CummPoison returns the area under the Poi(u)-distribution.                 }
{                                                                              }
{                                                                              }
function  CummPoisson(const X: Integer; const u: MFloat): MFloat;



{                                                                              }
{ TStatistic                                                                   }
{                                                                              }
{   Class that computes various descriptive statistics on a sample without     }
{   storing the sample values.                                                 }
{                                                                              }
{   To use, call one of the Add methods for every sample value. The values of  }
{   the descriptive statistics are available after every call to Add.          }
{                                                                              }
{   Statistics calculated:                                                     }
{   ---------------------                                                      }
{                                                                              }
{   Count is the number of sample values added.                                }
{                                                                              }
{   Sum is the sum of all sample values. SumOfSquares is the sum of the        }
{   squares of all sample values. Likewise SumOfCubes and SumOfQuads.          }
{                                                                              }
{   Min and Max are the minimum and maximum sample values. Range is the        }
{   difference between the maximum and minimum sample values.                  }
{                                                                              }
{   Mean (or average) is the sum of all data values divided by the number of   }
{   elements in the sample.                                                    }
{                                                                              }
{   Variance is a measure of the spread of a distribution about its mean and   }
{   is defined by var(X) = E([X - E(X)]2). The variance is expressed in the    }
{   squared unit of measurement of X.                                          }
{                                                                              }
{   Standard deviation is the square root of the variance and like variance    }
{   is a measure of variability or dispersion of a sample. Standard deviation  }
{   is expressed in the same unit of measurement as the sample values.         }
{   If a distribution's standard deviation is greater than its mean, the mean  }
{   is inadequate as a representative measure of central tendency. For         }
{   normally distributed data values, approximately 68% of the distribution    }
{   falls within ± 1 standard deviation of the mean and 95% of the             }
{   distribution falls within ± 2 standard deviations of the mean.             }
{                                                                              }
{   M1, M2, M3 and M4 are the first four central moments (moments about the    }
{   mean). The second moment about the mean is equal to the variance.          }
{                                                                              }
{   Skewness is the degree of asymmetry about a central value of a             }
{   distribution. A distribution with many small values and few large values   }
{   is positively skewed (right tail), the opposite (left tail) is negatively  }
{   skewed.                                                                    }
{                                                                              }
{   Kurtosis is the degree of peakedness of a distribution, defined as a       }
{   normalized form of the fourth central moment of a distribution. Kurtosis   }
{   is based on the size of a distribution's tails. Distributions with         }
{   relatively large tails are called "leptokurtic"; those with small tails    }
{   are called "platykurtic." A distribution with the same kurtosis as the     }
{   normal distribution is called "mesokurtic."  The kurtosis of a normal      }
{   distribution is 0.                                                         }
{                                                                              }

(*type
  TStatistic = class
  protected
    FCount        : Integer;
    FMin          : MFloat;
    FMax          : MFloat;
    FSum          : MFloat;
    FSumOfSquares : MFloat;
    FSumOfCubes   : MFloat;
    FSumOfQuads   : MFloat;

  public
    procedure Assign(const S: TStatistic);
    function  Duplicate: TStatistic;
    procedure Clear;
    function  IsEqual(const S: TStatistic): Boolean;

    procedure Add(const V: MFloat); overload;
    procedure Add(const V: Array of MFloat); overload;
    procedure Add(const V: TStatistic); overload;
    procedure AddNegated(const V: TStatistic);
    procedure Negate;

    property  Count: Integer read FCount;
    property  Min: MFloat read FMin;
    property  Max: MFloat read FMax;
    property  Sum: MFloat read FSum;
    property  SumOfSquares: MFloat read FSumOfSquares;
    property  SumOfCubes: MFloat read FSumOfCubes;
    property  SumOfQuads: MFloat read FSumOfQuads;

    function  Range: MFloat;
    function  Mean: MFloat;
    function  PopulationVariance: MFloat;
    function  PopulationStdDev: MFloat;
    function  Variance: MFloat;
    function  StdDev: MFloat;

    function  M1: MFloat;
    function  M2: MFloat;
    function  M3: MFloat;
    function  M4: MFloat;
    function  Skew: MFloat;
    function  Kurtosis: MFloat;

    function  GetAsString: String;
  end;
  EStatistic = class(EStatistics);
  EStatisticNoSample = class(EStatistic);
  EStatisticDivisionByZero = class(EStatistic);
     *)


{                                                                              }
{ Tests                                                                        }
{                                                                              }
{$IFDEF MATHS_TEST}
procedure Test;
{$ENDIF}



implementation

(*uses
  { System }
  Math,

  { Fundamentals }
  flcFloats;      *)



{                                                                              }
{ Binomial distribution                                                        }
{                                                                              }
function BinomialCoeff(N, R: Integer): MFloat;
var I, K : Integer;
begin
  if (N <= 0) or (R > N) then
    Xraise (EStatisticsInvalidArgument.Create('BinomialCoeff: Invalid argument'));
  if N > 1547 then
    Xraise (EStatisticsOverflow.Create('BinomialCoeff: Overflow'));
  Result := 1.0;
  if (R = 0) or (R = N) then
    exit;
  if R > N div 2 then
   R := N - R;
  K := 2;
  For I := N - R + 1 to N do begin
      Result := Result * I;
      if K <= R then begin
          Result := Result / K;
          Inc(K);
        end;
    end;
  Result := Int(Result + 0.5);
end;



{                                                                              }
{ gamma function, incomplete, series evaluation                                }
{ The gamma functions were translated from 'Numerical Recipes'.                }
{ 
(*                                                                             }
//const
 // {$IFDEF MFloatIsExtended}
 // CompareDelta = ExtendedCompareDelta;
 // {$ELSE}
 // {$IFDEF MFloatIsDouble}
  const CompareDelta = DoubleCompareDelta;
 // {$ENDIF}
 // {$ENDIF}   *)
 
 const itmax = 100;
      eps   = 3.0e-7;

procedure gser(const a, x: MFloat; var gamser, gln: MFloat);

var n : Integer;
    sum, del, ap : MFloat;
begin
  gln := GammaLn(a);
  if FloatZero(x, CompareDelta) then begin
      GamSer := 0.0;
      exit;
    end;
  if X < 0.0 then
    Xraise (EStatisticsInvalidArgument.Create('Gamma: GSER: Invalid argument'));
  ap := a;
  sum := 1.0 / a;
  del := sum;
  for n := 1 to itmax do begin
      ap := ap + 1.0;
      del := del * x / ap;
      sum := sum + del;
      if abs(del) < abs(sum) * eps then begin
          GamSer := sum * exp(-x + a * ln(x) - gln);
          exit;
        end;
    end;
  raise EStatisticsOverflow.Create('Gamma: GSER: Overflow: ' +
      'Argument a is too large or itmax is too small');
end;

{ gamma function, incomplete, continued fraction evaluation                    }
procedure gcf(const a, x: MFloat; var gammcf, gln: MFloat);
//const itmax = 100;
  //    eps   = 3.0e-7;
var n : integer;
    agold, g, fac, b1, b0, anf, ana, an, a1, a0 : MFloat;
begin
  gln := GammaLn(a);
  agold := 0.0;
  g := 0.0;
  a0 := 1.0;
  a1 := x;
  b0 := 0.0;
  b1 := 1.0;
  fac := 1.0;
  For n := 1 to itmax do begin
      an := 1.0 * n;
      ana := an - a;
      a0 := (a1 + a0 * ana) * fac;
      b0 := (b1 + b0 * ana) * fac;
      anf := an * fac;
      a1 := x * a0 + anf * a1;
      b1 := x * b0 + anf * b1;
      if not FloatZero(a1, CompareDelta) then begin
          fac := 1.0 / a1;
          g := b1 * fac;
          if abs((g - gold) / g) < eps then
            break;
          agold := g;
        end;
    end;
  Gammcf := exp(-x + a * ln(x) - gln) * g;
end;

{ GAMMP  gamma function, incomplete                                            }
function GammP(const a,x: MFloat): MFloat;
var gammcf, gln : MFloat;
begin
  if (x < 0.0) or (a <= 0.0) then
    Xraise (EStatisticsInvalidArgument.Create('Gamma: GAMMP: Invalid argument'));
  if x < a + 1.0 then begin
      gser(a, x, gammcf, gln);
      result := gammcf
    end
  else
    begin
      gcf(a, x, gammcf, gln);
      result := 1.0 - gammcf
    end;
end;

{ GAMMQ  gamma function, incomplete, complementary                             }
function gammq(const a, x: MFloat): MFloat;
var gamser, gln : MFloat;
begin
  if (x < 0.0) or (a <= 0.0) then
    Xraise( EStatisticsInvalidArgument.Create('Gamma: GAMMQ: Invalid argument'));
  if x < a + 1.0 then begin
      gser(a, x, gamser, gln);
      Result := 1.0 - gamser;
    end
  else begin
      gcf(a, x, gamser, gln);
      Result := gamser
    end;
end;

{ error function                                                               }
function erf(x: MFloat): MFloat;
begin
   if x < 0.0 then
     result := -gammp(0.5, sqr(x))
   else
     result := gammp(0.5, sqr(x))
end;

{ error function, complementary                                                }
function erfc(const x: MFloat): MFloat;
begin
  if x < 0.0 then
    result := 1.0 + gammp(0.5, sqr(x))
  else
    result := gammq(0.5, sqr(x));
end;



{                                                                              }
{ BETACF  beta function, incomplete, continued fraction evaluation             }
{ The beta functions were translated from 'Numerical Recipes'.                 }
{                                                                              }
function betacf(const a, b, x: MFloat): MFloat;
//const itmax = 100;
  //    eps   = 3.0e-7;
var tem, qap, qam, qab, em, d : MFloat;
    bz, bpp, bp, bm, az, app  : MFloat;
    am, aold, ap              : MFloat;
    m                         : Integer;
begin
  am := 1.0;
  bm := 1.0;
  az := 1.0;
  qab := a + b;
  qap := a + 1.0;
  qam := a - 1.0;
  bz := 1.0 - qab * x / qap;
  For m := 1 to itmax do begin
      em := m;
      tem := em + em;
      d := em * (b - m) * x / ((qam + tem) * (a + tem));
      ap := az + d * am;
      bp := bz + d * bm;
      d := -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
      app := ap + d * az;
      bpp := bp + d * bz;
      aold := az;
      am := ap / bpp;
      bm := bp / bpp;
      az := app / bpp;
      bz := 1.0;
      if abs(az - aold) < eps * abs(az) then begin
          Result := az;
          exit;
        end;
    end;
  raise EStatisticsOverflow.Create('Beta: BETACF: ' +
      'Argument a or b is too big or itmax is too small');
end;

{ BETAI  beta function, incomplete                                             }
function betai(const a, b, x: MFloat): MFloat;
var bt : MFloat;
begin
  if (x < 0.0) or (x > 1.0) then
    xraise (EStatisticsInvalidArgument.Create('Beta: BETAI: Invalid argument'));
    //SIRegister_cFundamentUtils(CL: TPSPascalCompiler);
  if FloatZero(x, CompareDelta) or FloatOne(x, CompareDelta) then
    bt := 0.0
  else
    bt := exp(GammaLn(a + b) - GammaLn(a) - GammaLn(b) + a * ln(x) + b * ln(1.0 - x));
  if x < (a + 1.0) / (a + b + 2.0) then
    Result := bt * betacf(a, b, x) / a
  else
    Result := 1.0 - bt * betacf(b, a, 1.0 - x) / b;
end;



{                                                                              }
{ Normal distribution                                                          }
{                                                                              }
function CummNormal(const u, s, X: MFloat): MFloat;
begin
  result := erfc(((X - u) / s) / Sqrt2) / 2.0;
end;

function CummNormal01(const X: MFloat): MFloat;
begin
  result := erfc(X / Sqrt2) / 2.0;
end;

{                                                                              }
{ Polynomial                                                                   }
{                                                                              }
function PolyEval(const X: MFloat; const Coef: array of MFloat; const N: Integer): MFloat;
var P : MFloat;
    L : Integer;
begin
  if Length(Coef) <> N + 1 then
    xraise (EInvalidArgument.Create('PolyEval: Invalid number of coefficients'));
  P := 1.0;
  L := N;
  Result := 0.0;
  while L >= 0 do begin
      Result := Result + Coef[L] * P;
      P := P * X;
      Dec(L);
    end;
end;

{                                                                   }
{ Returns the argument, x, for which the area under the             }
{ Gaussian probability density function (integrated from            }
{ minus infinity to x) is equal to y.                               }
{                                                                   }
{ For small arguments 0 < y < exp(-2), the program computes         }
{ z = sqrt( -2.0 * log(y) );  then the approximation is             }
{ x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).                        }
{ There are two rational functions P/Q, one for 0 < y < exp(-32)    }
{ and the other for y up to exp(-2).  For larger arguments,         }
{ w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).        }
{                                                                   }
{ Algorithm translated into Delphi by David Butler from Cephes C    }
{ library with persmission of Stephen L. Moshier                    }
{ <moshier@na-net.ornl.gov>.                                        }
{                                                                   }

const ExpM2        = 0.135335283236612691893999494972484; { e^-2   }

function InvCummNormal01(Y0: MFloat): MFloat;

var P0 : array[0..4] of MFloat;
    Q0 : array[0..8] of MFloat;
    P1 : array[0..8] of MFloat;
    Q1 : array[0..8] of MFloat;
    P2 : array[0..8] of MFloat;
    Q2 : array[0..8] of MFloat;
    
(*const P0 : array[0..4] of MFloat = (
           -5.99633501014107895267e1,
            9.80010754185999661536e1,
           -5.66762857469070293439e1,
            1.39312609387279679503e1,
           -1.23916583867381258016e0);
      Q0 : array[0..8] of MFloat = (
            1.00000000000000000000e0,
            1.95448858338141759834e0,
            4.67627912898881538453e0,
            8.63602421390890590575e1,
           -2.25462687854119370527e2,
            2.00260212380060660359e2,
           -8.20372256168333339912e1,
            1.59056225126211695515e1,
           -1.18331621121330003142e0);
      P1 : array[0..8] of MFloat = (
            4.05544892305962419923e0,
            3.15251094599893866154e1,
            5.71628192246421288162e1,
            4.40805073893200834700e1,
            1.46849561928858024014e1,
            2.18663306850790267539e0,
           -1.40256079171354495875e-1,
           -3.50424626827848203418e-2,
           -8.57456785154685413611e-4);
      Q1 : array[0..8] of MFloat = (
            1.00000000000000000000e0,
            1.57799883256466749731e1,
            4.53907635128879210584e1,
            4.13172038254672030440e1,
            1.50425385692907503408e1,
            2.50464946208309415979e0,
           -1.42182922854787788574e-1,
           -3.80806407691578277194e-2,
           -9.33259480895457427372e-4);
      P2 : array[0..8] of MFloat = (
            3.23774891776946035970e0,
            6.91522889068984211695e0,
            3.93881025292474443415e0,
            1.33303460815807542389e0,
            2.01485389549179081538e-1,
            1.23716634817820021358e-2,
            3.01581553508235416007e-4,
            2.65806974686737550832e-6,
            6.23974539184983293730e-9);
      Q2 : array[0..8] of MFloat = (
            1.00000000000000000000e0,
            6.02427039364742014255e0,
            3.67983563856160859403e0,
            1.37702099489081330271e0,
            2.16236993594496635890e-1,
            1.34204006088543189037e-2,
            3.28014464682127739104e-4,
            2.89247864745380683936e-6,
            6.79019408009981274425e-9); *)
var X, Z, Y2, X0, X1 : MFloat;
    Code             : Boolean;
begin
  if (Y0 <= 0.0) or (Y0 >= 1.0) then
    EStatisticsInvalidArgument.Create('InvCummNormal01: Invalid argument');
  Code := True;
  if Y0 > 1.0 - ExpM2 then //fix
    begin
      Y0 := 1.0 - Y0;
      Code := False;
    end;
   if Y0 > ExpM2 then begin
      Y0 := Y0 - 0.5;
      Y2 := Y0 * Y0;
      X := Y0 + Y0 * (Y2 * PolyEval(Y2, P0, 4) / PolyEval(Y2, Q0, 8));
      result := X * Sqrt2Pi;
    end
  else begin
      X := Sqrt(-2.0 * Ln(Y0));
      X0 := X - Ln(X) / X;
      Z := 1.0 / X;
      if X < 8.0 then
        X1 := Z * PolyEval(Z, P1, 8) / PolyEval(Z, Q1, 8) else
        X1 := Z * PolyEval(Z, P2, 8) / PolyEval(Z, Q2, 8);
      X := X0 - X1;
      if Code then
        X := -X;
      result := X;
    end;
end;

function InvCummNormal(const u, s, Y0: MFloat): MFloat;
begin
  result := InvCummNormal01(Y0) * s + u;
end;



{                                                                              }
{ Chi-Squared distribution                                                     }
{                                                                              }
function CummChiSquare(const Chi, Df: MFloat): MFloat;
begin
  result := 1.0 - gammq(0.5 * Df, 0.5 * Chi);
end;



{                                                                              }
{ F distribution                                                               }
{                                                                              }
function CumF(const f, Df1, Df2: MFloat): MFloat;
begin
  if F <= 0.0 then
    xraise (EStatisticsInvalidArgument.Create('CumF: Invalid argument'));
  result := 1.0 - (betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f))
       + (1.0 - betai(0.5 * df1, 0.5 * df2, df1 / (df1 + df2 / f)))) / 2.0;
end;



{                                                                              }
{ Poisson distribution                                                         }
{                                                                              }
function CummPoisson(const X: Integer; const u: MFloat): MFloat;
begin
  result := GammQ(X + 1, u);
end;

//////////////////// from flcMaths.pas /////////////////////////////////////////

procedure SetFPURoundingTruncate;
begin
  {$IFDEF DELPHI6_UP}
  SetRoundMode(rmTruncate);
  {$ENDIF}
end;

procedure PolarToRectangular(const R, Theta: MFloat; var X, Y: MFloat);
var S, C : MFloat;
begin
  glSinCos11(Theta, S, C);
  X := R * C;
  Y := R * S;
end;

procedure RectangularToPolar(const X, Y: MFloat; var R, Theta: MFloat);
begin
  if FloatZero(X, CompareDelta) then
    if FloatZero(Y, CompareDelta) then
      R := 0.0
    else
      if Y > 0.0 then
        R := Y
      else
        R := -Y
  else
    R := Sqrt(Sqr(X) + Sqr(Y));
  Theta := ArcTan2(Y, X);
end;

function Distance(const X1, Y1, X2, Y2: MFloat): MFloat;
begin
  Result := Sqrt(Sqr(X1 - X2) + Sqr(Y1 - Y2));
end;

function Sign(const R: Integer): Integer;
begin
  if R > 0 then
    Result := 1 else
    if R < 0 then
      Result := -1 else
      Result := 0;
end;
//{$ENDIF}

function Signint64(const R: Int64): Integer;
begin
  if R > 0 then
    Result := 1 else
    if R < 0 then
      Result := -1 else
      Result := 0;
end;

function SignSingle(const R: Single): Integer;
begin
  if R > 0.0 then
    Result := 1 else
    if R < 0.0 then
      Result := -1 else
      Result := 0;
end;

function SignDouble(const R: Double): Integer;
begin
  if R > 0.0 then
    Result := 1 else
    if R < 0.0 then
      Result := -1 else
      Result := 0;
end;

function SignExtended(const R: Extended): Integer;
begin
  if R > 0.0 then
    Result := 1 else
    if R < 0.0 then
      Result := -1 else
      Result := 0;
end;

function FloatMod(const A, B: MFloat): MFloat;
begin
  Result := A - Floor(A / B) * B;
end;

{                                                                              }
{ Trigonometric functions                                                      }
{                                                                              }
function ATan360(const X, Y: MFloat): MFloat;
var Angle: MFloat;
begin
  if FloatZero(X, CompareDelta) then
    Angle := PiOn2
  else
    Angle := ArcTan(Y / X);
  Angle := Angle * DegPerRad;
  if (X <= 0.0) and (Y < 0.0) then
    Angle := Angle - 180.0;
  if (X < 0.0) and (Y > 0.0) then
    Angle := Angle + 180.0;
  if Angle < 0.0 then
    Angle := Angle + 360.0;
  result := Angle;
end;

function InverseTangentDeg(const X, Y: MFloat): MFloat;
{ 0 <= Result <= 360 }
var Angle : MFloat;
begin
  if FloatZero(X, CompareDelta) then
    Angle := PiOn2
  else
    Angle := ArcTan (Y / X);
  Angle := Angle * 180.0 / Pi;
  if (X <= 0.0) and (Y < 0.0) then
    Angle := Angle - 180.0
  else
  if (X < 0.0) and (Y > 0.0) then
    Angle := Angle + 180.0;
  if Angle < 0.0 then
    Angle := Angle + 360.0;
  result := Angle;
end;

const Pi2          = 6.283185307179586476925286766559006; { Pi * 2                 }

function InverseTangentRad(const X, Y: MFloat): MFloat;
{ 0 <= result <= 2pi }
var Angle : MFloat;
begin
  if FloatZero(X, CompareDelta) then
    Angle := PiOn2
  else
    Angle := ArcTan(Y / X);
  if (X <= 0.0) and (Y < 0) then
    Angle := Angle - Pi;
  if (X < 0.0) and (Y > 0) then
    Angle := Angle + Pi;
  If Angle < 0 then
    Angle := Angle + Pi2;
  result := Angle;
end;

function InverseSinDeg(const Y, R: MFloat): MFloat;
{ -90 <= result <= 90 }
var X : MFloat;
begin
  X := Sqrt(Sqr(R) - Sqr(Y));
  Result := InverseTangentDeg(X, Y);
  If Result > 90.0 then
    Result := Result - 360.0;
end;

function InverseSinRad(const Y, R: MFloat): MFloat;
{ -90 <= result <= 90 }
var X : MFloat;
begin
  X := Sqrt(Sqr(R) - Sqr(Y));
  Result := InverseTangentRad(X, Y);
  if Result > 90.0 then
    Result := Result - 360.0;
end;

function InverseCosDeg(const X, R: MFloat): MFloat;
{ -90 <= result <= 90 }
var Y : MFloat;
begin
  Y := Sqrt(Sqr(R) - Sqr(X));
  Result := InverseTangentDeg(X, Y);
  if Result > 90.0 then
    Result := Result - 360.0;
end;

function InverseCosRad(const X, R: MFloat): MFloat;
{ -90 <= result <= 90 }
var Y : MFloat;
begin
  Y := Sqrt(Sqr(R) - Sqr(X));
  Result := InverseTangentRad(X, Y);
  if Result > 90.0 then
    Result := Result - 360.0;
end;

procedure RealToDMS(const X: MFloat; var Degs, Mins, Secs: MFloat);
var Y : MFloat;
begin
  Degs := Int(X);
  Y := Frac(X) * 60.0;
  Mins := Int(Y);
  Secs := Frac(Y) * 60.0;
end;

function DMSToReal(const Degs, Mins, Secs: MFloat): MFloat;
begin
  Result := Degs + Mins / 60.0 + Secs / 3600.0;
end;

function CanonicalForm(const Theta: MFloat): MFloat;                        {-PI < theta <= PI}
begin
  if Abs(Theta) > Pi then
     Result := Round(Theta / Pi2) * Pi2
  else
     Result := Theta;
end;

var
  WordPrimesInit : Boolean; //= False;
  WordPrimesSet  : array of Word32; // 8192 bytes when initialized, 1 bit for
                                      // each number between 0-65535.


{ IsPrimeFactor                                                                }
function IsPrimeFactor(const N, F: Int64): Boolean;
begin
  Result := (N mod F = 0) and IsPrime(F);
end;

{ PrimeFactors                                                                 }
function PrimeFactors(const N: Int64): TInt64DynArray; //'Int64Array;
var V : Int64;
    E : MFloat;
    M : Word32;
    I : Word32;
begin
  // Initialize
  Result := []; //nil;
  V := N;
  if V < 0 then
    V := -V;
  // 0 and 1 has no prime factors
  if V <= 1 then
    exit;
  // Calculate maximum factor
  E := V;
  M := Ceil(Sqrt(E));
  // Find prime factors
  for I := 2 to M do
    if IsPrime(I) and (V mod I = 0) then
      begin
        // I is a prime factor
        //DynArrayAppend(Result, I);
        repeat
          V := V div I;
          if V = 1 then
            // No more factors
            exit;
        until V mod I <> 0;
      end;
  // Check for remaining prime factor
  if IsPrime(V) then
    //DynArrayAppend(Result, V);
end;

{ Find the GCD using Euclid's algorithm                                        }
function GCD(const N1, N2: Integer): Integer;
var X, Y, J : Integer;
begin
  X := N1;
  Y := N2;
  if X < Y then
    synSwapint(X, Y);
  while (X <> 1) and (X <> 0) and (Y <> 1) and (Y <> 0) do begin
      J := (X - Y) mod Y;
      if J = 0 then
        begin
          Result := AbsInt(Y);
          exit;
        end;
      X := Y;
      Y := J;
    end;
  Result := 1;
end;

{ Find the modular inverse of A modulo N using Euclid's algorithm.             }
function InvMod(const A, N: Integer): Integer;
var g0, g1, v0, v1, y, z : Integer;
begin
  if N < 2 then
    xraise (EInvalidArgument.CreateFmt('InvMod: n=%d < 2', [N]));
  if GCD (A, N) <> 1 then
    xraise (EInvalidArgument.CreateFmt('InvMod: GCD (a=%d, n=%d) <> 1', [A, N])) ;
  g0 := N;  g1 := A;
  v0 := 0;  v1 := 1;
  while g1 <> 0 do begin
      y := g0 div g1;
      z := g1;
      g1 := g0 - y * g1;
      g0 := z;
      z := v1;
      v1 := v0 - y * v1;
      v0 := z;
    end;
  if v0 > 0 then
    Result := v0 else
    Result := v0 + N;
end;


{ Calculates x = a^z mod n                                                     }
function ExpMod(A, Z: Integer; const N: Integer): Integer;
var Signed : Boolean;
begin
  Signed := Z < 0;
  if Signed then
    Z := -Z;
  Result := 1;
  while Z <> 0 do
    begin
      while not Odd(Z) do begin
          Z := Z shr 1;
          A := (A * Int64(A)) mod N;
        end;
      Dec (Z);
      Result := (Result * Int64(A)) mod N;
    end;
  if Signed then
    Result := InvMod(Result, N);
end;

{                                                                              }
{ Numerical solvers                                                            }
{                                                                            }

type
  fx = function (const x: Float): Float;                                                                         
function SecantSolver(const f: fx; const y, Guess1, Guess2: MFloat): MFloat;
var xn, xnm1, xnp1, fxn, fxnm1 : MFloat;
begin
  xnm1 := Guess1;
  xn := Guess2;
  fxnm1 := f(xnm1) - y;
  repeat
    fxn := f(xn) - y;
    xnp1 := xn - fxn * (xn - xnm1) / (fxn - fxnm1);
    fxnm1 := fxn;
    xnm1 := xn;
    xn := xnp1;
  until (f(xn - 0.00000001) - y) * (f(xn + 0.00000001) - y) <= 0.0;
  Result := xn;
end;

function NewtonSolver(const f, df: fx; const y, Guess: MFloat): MFloat;
var xn, xnp1 : MFloat;
begin
  xnp1 := Guess;
  repeat
    xn := xnp1;
    xnp1 := xn - f(xn) / df(xn);
  until Abs(xnp1 - xn) < 0.000000000000001;
  Result := xn;
end;

const h = 1e-15;

function FirstDerivative(const f: fx; const x: MFloat): MFloat;
begin
  Result := (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
end;

function SecondDerivative(const f: fx; const x: MFloat): MFloat;
begin
  Result := (-f(x + 2 * h) + 16 * f(x + h) - 30 * f(x) +
             16 * f(x - h) - f(x - 2 * h)) / (12 * h * h);
end;

function ThirdDerivative(const f: fx; const x: MFloat): MFloat;
begin
  Result := (f(x + 2 * h) - 2 * f(x + h) + 2 * f(x - h) - f(x - 2 * h)) / (2 * h * h * h);
end;

function FourthDerivative(const f: fx; const x: MFloat): MFloat;
begin
  Result := (f(x + 2 * h) - 4 * f(x + h) + 6 * f(x) - 4 * f(x - h) + f(x - 2 * h)) / (h * h * h * h);
end;

function SimpsonIntegration(const f: fx; const a, b: MFloat; N: Integer): MFloat;
var hh : MFloat;
    I : Integer;
begin
  if N mod 2 = 1 then
   Inc(N); // N must be multiple of 2
  hh := (b - a) / N;
  Result := 0.0;
  for I := 1 to N - 1 do
    Result := Result + ((I mod 2) * 2 + 2) * f(a + (I - 0.5) * hh);
  Result := (Result + f(a) + f(b)) * hh / 3.0;
end;

{                                                                              }
{ Combinatorics                                                                }
{                                                                              }
function Combinations(const N, C: Integer): MFloat;
begin
  Result := Factorial(N) / (Factorial(C) * Factorial(N - C));
end;

function Permutations(const N, P: Integer): MFloat;
begin
  Result := Factorial(N) / Factorial(N - P);
end;



{                                                                              }
{ Fibonacci (N) = Fibonacci (N - 1) + Fibonacci (N - 2)                        }
{                                                                              }
function Fibonacci(const N: Integer): Int64;
var I      : Integer;
    f1, f2 : Int64;
begin
  if (N < 0) or (N > 92) then
    Xraise (EInvalidArgument.Create('Fibonacci: Invalid argument'));
  Result := 1;
  if N = 1 then
    exit;
  f1 := 0;     // fib(0) = 0
  f2 := 1;     // fib(1) = 1
  for I := 1 to N do begin
      Result := f1 + f2;
      f2 := f1;
      f1 := Result;
    end;
end;

procedure CalcFrequencyF(const FrequencyIndex: Integer;
    const RealIn, ImagIn: array of MFloat;
    var RealOut, ImagOut: MFloat);
var K, NumSamples                 : Integer;
    cos1, cos2, cos3, theta, beta : MFloat;
    sin1, sin2, sin3              : MFloat;
begin
  NumSamples := Length(RealIn);
  if NumSamples <> Length(ImagIn) then
    Xraise 
      (EInvalidArgument.Create('CalcFrequency: RealIn and ImagIn must be of equal length'));
  RealOut := 0.0;
  ImagOut := 0.0;
  theta := Pi2 * FrequencyIndex / NumSamples;
  sin1 := sin (-2 * theta);
  sin2 := sin (-theta);
  cos1 := cos (-2 * theta);
  cos2 := cos (-theta);
  beta := 2 * cos2;
  for K := 0 to NumSamples - 1 do begin
      sin3 := beta * sin2 - sin1;
      sin1 := sin2;
      sin2 := sin3;

      cos3 := beta * cos2 - cos1;
      cos1 := cos2;
      cos2 := cos3;

      RealOut := RealOut + RealIn[K] * cos3 - ImagIn[K] * sin3;
      ImagOut := ImagOut + ImagIn[K] * cos3 + RealIn[K] * sin3;
    end;
end;

//Procedure CalcFrequency(NumSamples,FrequencyIndex: Int;InArray: TCompVector;var FT : Complex)
//unit uPSI_JclStreams;




{                                                                              }
{ TStatistic                                                                   }
{                                                                              }
const
  StatisticFloatDelta = CompareDelta;

(*
procedure TStatisticAssign(const S: TStatistic);
begin
  Assert(Assigned(S), 'Assigned(S)');
  FCount := S.FCount;
  FMin := S.FMin;
  FMax := S.FMax;
  FSum := S.FSum;
  FSumOfSquares := S.FSumOfSquares;
  FSumOfCubes := S.FSumOfCubes;
  FSumOfQuads := S.FSumOfQuads;
end;

function TStatistic.Duplicate: TStatistic;
begin
  Result := TStatistic.Create;
  Result.Assign(self);
end;

procedure TStatistic.Clear;
begin
  FCount := 0;
  FMin := 0.0;
  FMax := 0.0;
  FSum := 0.0;
  FSumOfSquares := 0.0;
  FSumOfCubes := 0.0;
  FSumOfQuads := 0.0;
end;

function TStatistic.IsEqual(const S: TStatistic): Boolean;
begin
  Result :=
      Assigned(S) and
      (FCount = S.FCount) and
      (FMin = S.FMin) and
      (FMax = S.FMax) and
      (FSum = S.FSum) and
      (FSumOfSquares = S.FSumOfSquares) and
      (FSumOfCubes = S.FSumOfCubes) and
      (FSumOfQuads = S.FSumOfQuads);
end;

procedure TStatistic.Add(const V: MFloat);
var A: MFloat;
begin
  Inc(FCount);
  if FCount = 1 then
    begin
      FMin := V;
      FMax := V;
    end
  else
    begin
      if V < FMin then
        FMin := V
      else
        if V > FMax then
          FMax := V;
    end;
  FSum := FSum + V;
  A := Sqr(V);
  FSumOfSquares := FSumOfSquares + A;
  A := A * V;
  FSumOfCubes := FSumOfCubes + A;
  A := A * V;
  FSumOfQuads := FSumOfQuads + A;
end;

procedure TStatistic.Add(const V: Array of MFloat);
var I : Integer;
begin
  For I := 0 to High(V) - 1 do
    Add(V[I]);
end;

// Add the sample values of V
procedure TStatistic.Add(const V: TStatistic);
begin
  if Assigned(V) and (V.FCount > 0) then
    begin
      if FCount = 0 then
        begin
          FMin := V.FMin;
          FMax := V.FMax;
        end
      else
        begin
          if V.FMin < FMin then
            FMin := V.FMin;
          if V.FMax > FMax then
            FMax := V.FMax;
        end;
      Inc(FCount, V.FCount);
      FSum := FSum + V.FSum;
      FSumOfSquares := FSumOfSquares + V.FSumOfSquares;
      FSumOfCubes := FSumOfCubes + V.FSumOfCubes;
      FSumOfQuads := FSumOfQuads + V.FSumOfQuads;
    end;
end;

// Add the negated sample values of V
procedure TStatistic.AddNegated(const V: TStatistic);
begin
  if Assigned(V) and (V.FCount > 0) then
    begin
      Inc(FCount, V.FCount);
      if -V.FMax < FMin then
        FMin := -V.FMax;
      if -V.FMin > FMax then
        FMax := -V.FMin;
      FSum := FSum - V.FSum;
      FSumOfSquares := FSumOfSquares + V.FSumOfSquares;
      FSumOfCubes := FSumOfCubes - V.FSumOfCubes;
      FSumOfQuads := FSumOfQuads + V.FSumOfQuads;
    end;
end;

// Negate all sample values
procedure TStatistic.Negate;
var T: MFloat;
begin
  if FCount > 0 then
    begin
      T := FMin;
      FMin := -FMax;
      FMax := -T;
      FSum := -FSum;
      FSumOfCubes := -FSumOfCubes;
      // FSumOfSquares and FSumOfQuads stay unchanged
    end;
end;

function TStatistic.Range: MFloat;
begin
  if FCount > 0 then
    Result := FMax - FMin
  else
    Result := 0.0;
end;

function TStatistic.Mean: MFloat;
begin
  if FCount = 0 then
    raise EStatisticNoSample.Create('No mean');
  Result := FSum / FCount
end;

function TStatistic.PopulationVariance: MFloat;
begin
  if FCount > 0 then
    Result := (FSumOfSquares - Sqr(FSum) / FCount) / FCount
  else
    Result := 0.0;
end;

function TStatistic.PopulationStdDev: MFloat;
begin
  Result := Sqrt(PopulationVariance);
end;

// Variance is an unbiased estimator of s^2 (as opposed to PopulationVariance
// which is biased)
function TStatistic.Variance: MFloat;
begin
  if FCount > 1 then
    Result := (FSumOfSquares - Sqr(FSum) / FCount) / (FCount - 1)
  else
    Result := 0.0;
end;

function TStatistic.StdDev: MFloat;
begin
  Result := Sqrt(Variance);
end;

function TStatistic.M1: MFloat;
begin
  Result := FSum / (FCount + 1);
end;

function TStatistic.M2: MFloat;
var NI, M1 : MFloat;
begin
  NI := 1.0 / (FCount + 1);
  M1 := FSum * NI;
  Result := FSumOfSquares * NI - Sqr(M1);
end;

function TStatistic.M3: MFloat;
var NI, M1 : MFloat;
begin
  NI := 1.0 / (FCount + 1);
  M1 := FSum * NI;
  Result := FSumOfCubes * NI
          - M1 * 3.0 * FSumOfSquares * NI
          + 2.0 * Sqr(M1) * M1;
end;

function TStatistic.M4: MFloat;
var NI, M1, M1Sqr : MFloat;
begin
  NI := 1.0 / (FCount + 1);
  M1 := FSum * NI;
  M1Sqr := Sqr(M1);
  Result := FSumOfQuads * NI
          - M1 * 4.0 * FSumOfCubes * NI
          + M1Sqr * 6.0 * FSumOfSquares * NI
          - 3.0 * Sqr(M1Sqr);
end;

function TStatistic.Skew: MFloat;
begin
  Result := M3 * Power(M2, -3/2);
end;

function TStatistic.Kurtosis: MFloat;
var M2Sqr : MFloat;
begin
  M2Sqr := Sqr(M2);
  if FloatZero(M2Sqr, StatisticFloatDelta) then
    raise EStatisticDivisionByZero.Create('Kurtosis: Division by zero');
  Result := M4 / M2Sqr;
end;

function TStatistic.GetAsString: String;
const
  NL = #13#10;
begin
  if Count > 0 then
    Result := 'n: ' + IntToStr(Count) +
              '  Sum: ' + FloatToStr(Sum) +
              '  Sum of squares: ' + FloatToStr(SumOfSquares) +
              NL +
              'Sum of cubes: ' + FloatToStr(SumOfCubes) +
              '  Sum of quads: ' + FloatToStr(SumOfQuads) +
              NL +
              'Min: ' + FloatToStr(Min) +
              '  Max: ' + FloatToStr(Max) +
              '  Range: ' + FloatToStr(Range) +
              NL +
              'Mean: ' + FloatToStr(Mean) +
              '  Variance: ' + FloatToStr(Variance) +
              '  Std Dev: ' + FloatToStr(StdDev) +
              NL +
              'M3: ' + FloatToStr(M3) +
              '  M4: ' + FloatToStr(M4) +
              NL +
              'Skew: ' + FloatToStr(Skew) +
              '  Kurtosis: ' + FloatToStr(Kurtosis) +
              NL
  else
    Result := 'n: 0';
end;

  *)

{                                                                              }
{ Tests                                                                        }
{                                                                              }
{$IFDEF MATHS_TEST}
//{$ASSERTIONS ON}
procedure Test;
var A, B : TStatistic;
begin
  A := TStatistic.Create;
  B := TStatistic.Create;

  Assert(A.Count = 0);
  Assert(A.Sum = 0.0);

  A.Add(1.0);
  A.Add(2.0);
  A.Add(3.0);

  Assert(A.Count = 3);
  Assert(A.Sum = 6.0);
  Assert(A.Min = 1.0);
  Assert(A.Max = 3.0);
  Assert(A.Mean = 2.0);

  B.Assign(A);
  Assert(B.Sum = 6.0);

  B.Add(A);
  Assert(B.Sum = 12.0);

  A.Clear;
  Assert(A.Count = 0);

  A.Add(4.0);
  A.Add(10.0);
  A.Add(1.0);

  Assert(A.Count = 3);
  Assert(A.Sum = 15.0);
  Assert(A.Min = 1.0);
  Assert(A.Max = 10.0);
  Assert(A.Mean = 5.0);
  Assert(A.GetAsString <> '');

  B.Free;
  A.Free;
end;
{$ENDIF}

function ReverseBits32(const Value: Word32; const BitCount: Integer): Word32;
var I, C : Integer;
    V : Word32;
begin
  V := Value;
  Result := 0;
  C := BitCount;
  if C > 32 then
    C := 32;
  for I := 0 to C - 1 do begin
      Result := (Result shl 1) or (V and 1);
      V := V shr 1;
    end;
end;

function RotateLeftBits16(const Value: Word; const Bits: Byte): Word;
var I, B : Integer;
    R : Word;
begin
  R := Value;
  if Bits >= 16 then
    B := Bits mod 16
  else
    B := Bits;
  for I := 1 to B do
    if R and $8000 = 0 then
      R := Word(R shl 1)
    else
      R := Word(R shl 1) or 1;
  Result := R;
end;

function RotateLeftBits32(const Value: Word32; const Bits: Byte): Word32;
var I, B : Integer;
    R : Word32;
begin
  R := Value;
  if Bits >= 32 then
    B := Bits mod 32
  else
    B := Bits;
  for I := 1 to B do
    if R and $80000000 = 0 then
      R := Word32(R shl 1)
    else
      R := Word32(R shl 1) or 1;
  Result := R;
end;

function RotateRightBits32(const Value: Word32; const Bits: Byte): Word32;
var I, B : Integer;
    R : Word32;
begin
  R := Value;
  if Bits >= 32 then
    B := Bits mod 32
  else
    B := Bits;
  for I := 1 to B do
    if R and 1 = 0 then
      R := Word32(R shr 1)
    else
      R := Word32(R shr 1) or $80000000;
  Result := R;
end;


{                                                                              }
{ ASCII functions                                                              }
{                                                                              }
function AsciiHexCharValueB(const C: Char): Integer;
begin
  case Ord(C) of
    Ord('0')..Ord('9') : Result := Ord(C) - Ord('0');
    Ord('A')..Ord('F') : Result := Ord(C) - Ord('A') + 10;
    Ord('a')..Ord('f') : Result := Ord(C) - Ord('a') + 10;
  else
    Result := -1;
  end;
end;

{                                                                              }
{ ASCII case conversion                                                        }
{     

                                                                         }
 const
  AsciiDigit0 = (#48);
  AsciiDigit9 = (#57);
  AsciiUpperA = (#65);
  AsciiUpperZ = (#90);
  AsciiLowerA = (#97);
  AsciiLowerZ = (#122);

//const
  //AsciiCaseDiff = Byte(AsciiLowerA) - Byte(AsciiUpperA);


var myc: TCompVector; //= array of Complex;
   var aFT : Complex;    //ast: TStringBuilder;
   
var myc2r, myc2i:  array of double; 
   var RealOut, ImagOut: MFloat;

begin

//function CummChiSquare(const Chi, Df: MFloat): MFloat;

 writeln(floattostr( CummChiSquare(73.5, 90.5)));
 
 Assert(not IsPrime(0),       'IsPrime(0)');
  Assert(not IsPrime(1),       'IsPrime(1)');
  Assert(IsPrime(2),           'IsPrime(2)');
  Assert(IsPrime(3),           'IsPrime(3)');
  //Assert2(IsPrime(-3),          'IsPrime(-3)');
  Assert(not IsPrime(4),       'IsPrime(4)');
  Assert(not IsPrime(-4),      'IsPrime(-4)');
  Assert(IsPrime(257),         'IsPrime(257)');
  Assert(not IsPrime($10002),  'IsPrime($10002)');
  Assert(IsPrime($10003),      'IsPrime($10003)');
  
  {type Complex = record //for complex numbers
    X, Y : Float;
  end;
  type TComplex', 'record Form : ComplexForm; X : Float; Y : Float; R :'
   + Float; Theta : Float; end;
  type TVector     = array of Float;
  TIntVector  = array of Int;
  TCompVector = array of Complex;}
  
  setlength(myc,2)
  
  myc[0].x:= 7.8   myc[1].x:= 17.8
  myc[0].y:= 6.7   myc[1].y:= 16.7
  
  
  CalcFrequency(2,0, myc, aFT)
  writeln(floattostr(aFT.X))
  writeln(floattostr(aFT.Y))
  
  //procedure CalcFrequencyF(const FrequencyIndex: Integer;
    //const RealIn, ImagIn: array of MFloat;
    //var RealOut, ImagOut: MFloat);
    
    setlength(myc2r,2)
    setlength(myc2i,2)
    myc2r[0]:= 7.8   myc2i[0]:= 17.8
    myc2r[1]:= 6.7   myc2i[1]:= 16.7
    CalcFrequencyF(0, myc2r, myc2i, RealOut, ImagOut);
    
    writeln(floattostr(RealOut))
    writeln(floattostr(imagOut))
    //const RealIn, ImagIn: array of MFloat;
    //var RealOut, ImagOut: MFloat);



End.

(*
{{                                                                              }
{ ByteChar                                                                     }
{   ByteChar is an one byte character.                                         }
{                                                                              }
{$IFDEF SupportAnsiChar}
type
  ByteChar = AnsiChar;
{$ELSE}
{$IFDEF SupportUTF8Char}
type
  ByteChar = UTF8Char;
{$ELSE}
type
  ByteChar = Byte;
  {$DEFINE ByteCharIsOrd}
{$ENDIF}
{$ENDIF}
type
  PByteChar = ^ByteChar;       *)

