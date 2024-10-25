with Symmetry_Operations; use Symmetry_Operations;

package Shmueli_Symbol_Parser is
   
   function Decode_Shmueli_Symbol (Symbol : in String)
                                  return Symmetry_Operator_Array;
   
end Shmueli_Symbol_Parser;
