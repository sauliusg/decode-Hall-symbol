with Symmetry_Operations; use Symmetry_Operations;

package Hall_Symbol_Parser is
   
   Debug_Print_Matrices : Boolean := False;
   
   function Decode_Hall_Symbol (Symbol : in String)
                               return Symmetry_Operator_Array;
   
end Hall_Symbol_Parser;
