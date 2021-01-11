{******************************************************************************}
{                                                                              }
{   Library:          Fundamentals 5.00                                        }
{   File name:        flcSysUtils.pas                                          }
{   File version:     5.02                                                     }
{   Description:      System utility functions                                 }
{                                                                              }
{   Copyright:        Copyright (c) 1999-2020, David J Butler                  }
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
{   2018/08/13  5.01  Initial version from other units.                        }
{   2019/07/29  5.02  FPC/Linux fixes.                                         }
{                                                                              }
{******************************************************************************}

{$INCLUDE ..\flcInclude.inc}

unit flcSysUtils;

interface

 type

  TfreqObj=class(TObject)
     ftemp,f,P,a,shape:INTEGER;
     StringRep:String;
     constructor Create(newf, newP, newA, newshape:integer);
     procedure makestringrep;
  end;
  
  
  type 
    MSArray = array[0..1] of tmemorystream;
  
  var 
  MS: MSArray;
    StreamInUse:integer;



function GetLastOSErrorCode: NativeInt;
function GetLastOSErrorMessage: String;
procedure TForm1msPlaySound(MS: MSArray; streaminuse: integer; aloop: boolean);
{                                                                              }
{ TBytes functions                                                             }
{                                                                              }
procedure BytesSetLengthAndZero(var V: TBytes; const NewLength: NativeInt);

procedure BytesInit(var V: TBytes; const R: Byte); //overload;
procedure BytesInit1(var V: TBytes; const S: String); //overload;

function  BytesAppend(var V: TBytes; const R: Byte): NativeInt; //overload;
function  BytesAppend1(var V: TBytes; const R: TBytes): NativeInt; //overload;
function  BytesAppend2(var V: TBytes; const R: array of Byte): NativeInt; //overload;
function  BytesAppend3(var V: TBytes; const R: String): NativeInt; //overload;

function  BytesCompare(const A, B: TBytes): Integer;

function  BytesEqual(const A, B: TBytes): Boolean;

 


implementation

uses
  {$IFDEF MSWINDOWS}
  Windows,
  {$ENDIF}

  {$IFDEF POSIX}
  {$IFDEF DELPHI}
  Posix.Errno,
  Posix.Unistd,
  {$ENDIF}
  {$ENDIF}

  {$IFDEF FREEPASCAL}
  {$IFDEF UNIX}
  BaseUnix,
  Unix,
  {$ENDIF}
  {$ENDIF}

  SysUtils;



{$IFDEF MSWINDOWS}
function GetLastOSErrorCode: NativeInt;
begin
  Result := NativeInt(Windows.GetLastError);
end;
{$ENDIF}

{$IFDEF DELPHI}
{$IFDEF POSIX}
function GetLastOSErrorCode: NativeInt;
begin
  Result := NativeInt(GetLastError);
end;
{$ENDIF}
{$ENDIF}

{$IFDEF FREEPASCAL}
{$IFDEF UNIX}
function GetLastOSErrorCode: NativeInt;
begin
  Result := NativeInt(GetLastOSError);
end;
{$ENDIF}
{$ENDIF}


  {************ TFreqObj.Create *********}
Constructor TFreqObj.Create;

begin
  inherited create;
  f:=newf;
  if f=0 then f:=1;
  ftemp:=f;
  p:=newp;
  a:=newA;
  shape:=newSHAPE;
  makestringrep;
end;

procedure TFreqObj.makestringrep;
begin
  stringrep:=format('%5d (F:%4d, P:%4d, A:%4d  %S)',[ftemp,f,p,a,SHAPENAMES[SHAPE]]);
end;

{********** msPlaySound ********}
 procedure TForm1msPlaySound(MS: MSArray; streaminuse: integer; aloop: boolean);
 var options:integer;
 begin
   options:=SND_MEMORY or SND_ASYNC;
   if aloop then  options:=options or SND_LOOP;
   PlaySound(MS[streaminuse].Memory, 0, options);
 end;

 {************** SaveSound **********}
Procedure TForm1savesound;
var
  ff:textfile;
  i:integer;
begin
  assignfile(ff,soundfilename);
  rewrite(ff);
  writeln(ff,'D',' ',duration.position,' ',unitsgrp.itemindex);
  for i:=0 to listbox1.items.count-1 do
  with listbox1, TFreqObj(items.objects[i]) do writeln(ff,
          integer(checked[i]),' ', f,' ',p,' ',a,' ',shape,' ',stringrep);
  closefile(ff);
  modified:=false;
end;

{                                                                              }
{ TBytes functions                                                             }
{                                                                              }
procedure BytesSetLengthAndZero(var V: TBytes; const NewLength: NativeInt);
var
  OldLen, NewLen : NativeInt;
begin
  NewLen := NewLength;
  if NewLen < 0 then
    NewLen := 0;
  OldLen := Length(V);
  if OldLen = NewLen then
    exit;
  SetLength(V, NewLen);
  if OldLen > NewLen then
    exit;
  ZeroMem(V[OldLen], NewLen - OldLen);
end;

procedure BytesInit(var V: TBytes; const R: Byte);
begin
  SetLength(V, 1);
  V[0] := R;
end;

procedure BytesInit1(var V: TBytes; const S: String);
var L, I : Integer;
begin
  L := Length(S);
  SetLength(V, L);
  for I := 0 to L - 1 do
    V[I] := Ord(S[I + 1]);
end;

function BytesAppend(var V: TBytes; const R: Byte): NativeInt;
begin
  Result := Length(V);
  SetLength(V, Result + 1);
  V[Result] := R;
end;

function BytesAppend1(var V: TBytes; const R: TBytes): NativeInt;
var
  L : NativeInt;
begin
  Result := Length(V);
  L := Length(R);
  if L > 0 then
    begin
      SetLength(V, Result + L);
      MoveMem(R[0], V[Result], L);
    end;
end;

function BytesAppend2(var V: TBytes; const R: array of Byte): NativeInt;
var
  L : NativeInt;
begin
  Result := Length(V);
  L := Length(R);
  if L > 0 then
    begin
      SetLength(V, Result + L);
      MoveMem(R[0], V[Result], L);
    end;
end;

function BytesAppend3(var V: TBytes; const R: String): NativeInt;
var
  L, I : NativeInt;
begin
  Result := Length(V);
  L := Length(R);
  if L > 0 then
    begin
      SetLength(V, Result + L);
      for I := 1 to L do
        V[Result] := Ord(R[I]);
    end;
end;

function BytesCompare(const A, B: TBytes): Integer;
var
  L, N : NativeInt;
begin
  L := Length(A);
  N := Length(B);
  if L < N then
    Result := -1
  else
  if L > N then
    Result := 1
  else
    Result := CompareMemB(Pointer(A)^, Pointer(B)^, L);
end;

function BytesEqual(const A, B: TBytes): Boolean;
var
  L, N : NativeInt;
begin
  L := Length(A);
  N := Length(B);
  if L <> N then
    Result := False
  else
    Result := EqualMem(Pointer(A)^, Pointer(B)^, L);
end;





resourcestring
  SSystemError = 'System error #%s';

{$IFDEF MSWINDOWS}
{$IFDEF StringIsUnicode}
function GetLastOSErrorMessage: String;
const MAX_ERRORMESSAGE_LENGTH = 256;
var Err: LongWord;
    Buf: array[0..MAX_ERRORMESSAGE_LENGTH - 1] of Word;
    Len: LongWord;
begin
  Err := Windows.GetLastError;
  FillChar(Buf, Sizeof(Buf), #0);
  Len := Windows.FormatMessageW(FORMAT_MESSAGE_FROM_SYSTEM, nil, Err, 0,
      @Buf, MAX_ERRORMESSAGE_LENGTH, nil);
  if Len = 0 then
    Result := Format(SSystemError, [IntToStr(Err)])
  else
    Result := StrPas(PWideChar(@Buf));
end;
{$ELSE}
function GetLastOSErrorMessage: String;
const MAX_ERRORMESSAGE_LENGTH = 256;
var Err: LongWord;
    Buf: array[0..MAX_ERRORMESSAGE_LENGTH - 1] of Byte;
    Len: LongWord;
begin
  Err := Windows.GetLastError;
  FillChar(Buf, Sizeof(Buf), #0);
  Len := Windows.FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM, nil, Err, 0,
      @Buf, MAX_ERRORMESSAGE_LENGTH, nil);
  if Len = 0 then
    Result := Format(SSystemError, [IntToStr(Err)])
  else
    Result := StrPas(PAnsiChar(@Buf));
end;
{$ENDIF}
{$ENDIF}

{$IFDEF DELPHI}
{$IFDEF POSIX}
function GetLastOSErrorMessage: String;
begin
  Result := SysErrorMessage(GetLastError);
end;
{$ENDIF}
{$ENDIF}

{$IFDEF FREEPASCAL}
{$IFDEF UNIX}
function GetLastOSErrorMessage: String;
begin
  Result := SysErrorMessage(GetLastOSError);
end;
{$ENDIF}
{$ENDIF}



end.

