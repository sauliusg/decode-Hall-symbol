with Symmetry_Operations; use Symmetry_Operations;

package HM_Symbol_Parser is
   
   Debug_Print_Matrices : Boolean := False;
   
   HM_SYMBOL_NOT_FOUND : exception;
   
   function Decode_Hermann_Mauguin_Symbol (Symbol : in String)
                                          return Symmetry_Operator_Array;
   
end HM_Symbol_Parser;
