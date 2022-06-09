with Text_IO;             use Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Ada.Command_Line;    use Ada.Command_Line;

procedure Decode_Hall is
   
   UNKNOWN_AXIS : exception;
   UNKNOWN_ROTATION : exception;
   UNKNOWN_TRANSLATION : exception;
   
   type Symop is array (1..4, 1..4) of Float;
   
   type Symop_Array is array (Positive range <>) of Symop;
   
   Zero_Matrix : constant Symop := (others => (others => 0.0));
   
   Unity_Matrix : constant Symop :=
     (
      (1 => 1.0, others => 0.0),
      (2 => 1.0, others => 0.0),
      (3 => 1.0, others => 0.0),
      (4 => 1.0, others => 0.0)
     );
   
   Ci_Matrix : constant Symop :=
     (
      (1 => -1.0, others => 0.0),
      (2 => -1.0, others => 0.0),
      (3 => -1.0, others => 0.0),
      (4 => -1.0, others => 0.0)
     );
   
   Inversion_Matrices : constant array (1..2) of Symop :=
     (Unity_Matrix, Ci_Matrix);
   
   type Crystallographic_Translation_Component is record
      Numerator : Integer range 0..6;
      Denominator : Integer range 1..6;
   end record;
      
   type Crystallographic_Translation is array (1..3)
     of Crystallographic_Translation_Component;
   
   A_Translation_Vector : constant Crystallographic_Translation :=
     ((0,1), (1,2), (1,2));
   
   B_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (0,1), (1,2));
   
   C_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (1,2), (0,1));
   
   I_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (1,2), (1,2));
   
   R_Translation_Vector_1 : constant Crystallographic_Translation :=
     ((1,3), (2,3), (2,3));
   
   R_Translation_Vector_2 : constant Crystallographic_Translation :=
     ((2,3), (1,3), (1,3));
   
   F_Translation_Vector_1 : constant Crystallographic_Translation :=
     A_Translation_Vector;
   
   F_Translation_Vector_2 : constant Crystallographic_Translation :=
     B_Translation_Vector;
   
   F_Translation_Vector_3 : constant Crystallographic_Translation :=
     C_Translation_Vector;
   
   
   Translation_a : constant Crystallographic_Translation :=
     ((1,2), (0,1), (0,1));
   Translation_b : constant Crystallographic_Translation :=
     ((0,1), (1,2), (0,1));
   Translation_c : constant Crystallographic_Translation :=
     ((0,1), (0,1), (1,2));
   Translation_n : constant Crystallographic_Translation :=
     ((1,2), (1,2), (1,2));
   Translation_u : constant Crystallographic_Translation :=
     ((1,4), (0,1), (0,1));
   Translation_v : constant Crystallographic_Translation :=
     ((0,1), (1,4), (0,1));
   Translation_w : constant Crystallographic_Translation :=
     ((0,1), (0,1), (1,4));
   Translation_d : constant Crystallographic_Translation :=
     ((1,4), (1,4), (1,4));
   
   Translations_3_1 : constant Crystallographic_Translation_Component := (1,3);
   Translations_3_2 : constant Crystallographic_Translation_Component := (2,3);
   
   Translations_4_1 : constant Crystallographic_Translation_Component := (1,4);
   Translations_4_3 : constant Crystallographic_Translation_Component := (3,4);
   
   Translations_6_1 : constant Crystallographic_Translation_Component := (1,6);
   Translations_6_2 : constant Crystallographic_Translation_Component := (2,6);
   Translations_6_4 : constant Crystallographic_Translation_Component := (4,6);
   Translations_6_5 : constant Crystallographic_Translation_Component := (5,6);
   
   procedure Put (S : Symop) is
   begin
      for J in Symop'Range(1) loop
         for K in Symop'Range(1) loop
            Put (" " & S (J,K)'Image);
         end loop;
         New_Line;
      end loop;
   end;
   
   procedure Put (S : Crystallographic_Translation) is
      Ratio : Float;
   begin
      for I in S'Range loop
         Ratio := Float (S (I).Numerator) / Float (S (I).Denominator); 
         Put (" " & Ratio'Image);
      end loop;
      New_Line;
   end;
   
   function To_Symop (T : Crystallographic_Translation) return Symop is
      S : Symop := Zero_Matrix;
   begin
      for I in T'Range loop
         S (I,4) := Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      return S;
   end;
   
   function To_Symop (T : Crystallographic_Translation_Component;
                      Axis : Positive) return Symop is
      S : Symop := Zero_Matrix;
   begin
      S (Axis,4) := Float (T.Numerator) / Float (T.Denominator);
      return S;
   end;
   
   Principal_Rotations : constant array (1..3, 1..4) of Symop :=
     (
      -- axis x (a)
      1 => (
            1 => ( -- twofold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0, -1.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            2 => ( --treefold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  1.0, -1.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            3 => ( -- fourfold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  1.0,  0.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            4 => ( --sixfold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  1.0, -1.0,  0.0),
                  (0.0,  1.0,  0.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 )
           ),
      -- axis y (b)
      2 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  (-1.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            3 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            4 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      -- axis z (c)
      3 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            3 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            4 => (
                  ( 1.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           )

     );
   
   Face_Diagonal_Rotations : constant array (1..3,1..2) of Symop :=
     (
      1 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      2 => (
            1 => (
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      3 => (
            1 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           )
     );
   
   Body_Diagonal_Rotation : constant Symop :=
     (
      (0.0, 0.0, 1.0, 0.0),
      (1.0, 0.0, 0.0, 0.0),
      (0.0, 1.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 1.0)
     );
   
   procedure Add (M1 : in out Symop; M2 : in Symop) is
   begin
      for I in M1'Range(1) loop
         for J in M1'Range(2) loop
            M1 (I,J) := M1 (I,J) + M2(I,J);
         end loop;
      end loop;
   end;
   
   function "+" (M1, M2 : Symop) return Symop is
      M : Symop;
   begin
      for I in M1'Range(1) loop
         for J in M1'Range(2) loop
            M (I,J) := M1 (I,J) + M2(I,J);
         end loop;
      end loop;
      return M;
   end;
   
   function "*" (M1, M2 : Symop) return Symop is
      M : Symop;
   begin
      pragma Assert (M1'Last(2) = M2'Last(1));
      
      for I in M1'Range(1) loop
         for J in M2'Range(2) loop
            M (I,J) := 0.0;
            for K in M2'Range(1) loop
               M (I,J) := M (I,J) + M1 (I,K) * M2(K,J);
            end loop;
         end loop;
      end loop;
      return M;
   end;
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
   procedure Get_Hall_Symbol_Inversions
     (
      Symbol : in String;
      Pos : in out Integer;
      N_Inversions : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      if Symbol (Pos) = '-' then
         N_Inversions := 2;
         Pos := Pos + 1;
      else
         N_Inversions := 1;
      end if;
   end;
   
   procedure Get_Hall_Symbol_Centerings
     (
      Symbol : in String;
      Pos : in out Positive;
      Centering : out Symop_Array;
      N_Centering : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      Centering (1) := Zero_Matrix;
      case Symbol (Pos) is
         when 'P' =>
           N_Centering := 1;
         when 'A' =>
           Centering (2) := To_Symop (A_Translation_Vector);
           N_Centering := 2;
         when 'B' =>
           Centering (2) := To_Symop (B_Translation_Vector);
           N_Centering := 2;
         when 'C' =>
           Centering (2) := To_Symop (C_Translation_Vector);
           N_Centering := 2;
         when 'I' =>
           Centering (2) := To_Symop (I_Translation_Vector);
           N_Centering := 2;
         when 'F' =>
           Centering (2) := To_Symop (F_Translation_Vector_1);
           Centering (3) := To_Symop (F_Translation_Vector_2);
           Centering (4) := To_Symop (F_Translation_Vector_3);
           N_Centering := 4;
         when 'R' =>
           Centering (2) := To_Symop (R_Translation_Vector_1);
           Centering (3) := To_Symop (R_Translation_Vector_2);
           N_Centering := 3;
         when others =>
           Centering (1) := Unity_Matrix;
           N_Centering := 1;           
      end case;
      Pos := Pos + 1;
   end;
   
   procedure Get_Hall_Symbol_Rotations
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotations : out Symop_Array;
      N_Rotations : in out Natural;
      Preceeding_Axis : in out Natural;
      Axis_Number : in Positive
     )
   is
      
      procedure Get_Inversion_Character (Inversion : out Character) is
      begin
         if Pos <= Symbol'Last and then
           Symbol (Pos) = '-' then
            Inversion := Symbol (Pos);
            Pos := Pos + 1;
         else
            Inversion := ' ';
         end if;
      end;
      
      procedure Get_Rotation_Character (Rotation : out Character) is
      begin
         if Pos <= Symbol'Last then
            pragma Assert (Symbol (Pos) in '2'..'6');
            Rotation := Symbol (Pos);
            Pos := Pos + 1;
         else
            Rotation := ' ';
         end if;
      end;
      
      procedure Get_Axis_Character (Axis : out Character;
                                    Axis_Number : Positive) is
      begin
         if Pos <= Symbol'Last and then
           (
            Symbol (Pos) in 'x'..'z' or else
              Symbol (Pos) = ''' or else
              Symbol (Pos) = '*' or else
              Symbol (Pos) = '"'
           )
         then
            Axis := Symbol (Pos);
            Pos := Pos + 1;
         else
            case Axis_Number is
               when 1 => Axis := 'x';
               when 2 => Axis := 'y';
               when 3 => Axis := 'z';
               when others =>
                  raise UNKNOWN_AXIS with "axis number" & Axis_Number'Image;
            end case;
         end if;
      end;
      
      procedure Get_Translation_Characters (Translations : out String) is
         I : Positive := 1;
      begin
         while Pos <= Symbol'Last and then
           (
            Symbol (Pos) in 'a' .. 'd' or else
              Symbol (Pos) in 'u' .. 'w' or else
              Symbol (Pos) in '1' .. '5' or else
              Symbol (Pos) = 'n' 
           ) loop
            Translations (I) := Symbol (Pos);
            Pos := Pos + 1;
            I := I + 1;
         end loop;
      end;
      
      function Rotation_Axis_Index (Rotation_Character : Character)
                                   return Positive
      is
      begin
         case Rotation_Character is
            when '2' => return 1;
            when '3' => return 2;
            when '4' => return 3;
            when '6' => return 6;
            when others => 
               raise UNKNOWN_ROTATION
                 with "rotation " & Rotation_Character'Image;
         end case;
      end;
      
      procedure Construct_Rotation_Matrix
        (
         Matrix : out Symop;
         Inversion : in Character;
         Axis : in Character;
         Rotation : in Character;
         Translations : String;
         Preceeding_Axis : in out Natural
        )
      is
         Axis_Number : Integer range 0..3 := 0;
      begin
         case Axis is 
            when 'x' =>
               Matrix :=
                 Principal_Rotations (1, Rotation_Axis_Index (Rotation));
               Preceeding_Axis := 1;
               Axis_Number := 1;
            when 'y' =>
               Matrix :=
                 Principal_Rotations (2, Rotation_Axis_Index (Rotation));
               Preceeding_Axis := 2;
               Axis_Number := 2;
            when 'z' =>
               Matrix :=
                 Principal_Rotations (3, Rotation_Axis_Index (Rotation));
               Axis_Number := 3;
            when ''' =>
               Matrix :=
                 Face_Diagonal_Rotations (Preceeding_Axis, 1);
               Axis_Number := 0;
            when '"' =>
               Matrix :=
                 Face_Diagonal_Rotations (Preceeding_Axis, 2);
               Axis_Number := 0;
            when '*' =>
               Matrix :=
                 Body_Diagonal_Rotation;
               Axis_Number := 0;
            when others =>
               raise UNKNOWN_AXIS with "axis character " & Axis'Image;
         end case;
         
         Preceeding_Axis := Axis_Number;
         
         if Inversion = '-' then
            Matrix := Ci_Matrix * Matrix;
         end if;
         
         for Tr of Translations loop
            case Tr is
               when ' ' => null;
               when 'a' => Add (Matrix, To_Symop (Translation_a));
               when 'b' => Add (Matrix, To_Symop (Translation_b));
               when 'c' => Add (Matrix, To_Symop (Translation_c));
               when 'd' => Add (Matrix, To_Symop (Translation_d));
               when 'u' => Add (Matrix, To_Symop (Translation_u));
               when 'v' => Add (Matrix, To_Symop (Translation_v));
               when 'w' => Add (Matrix, To_Symop (Translation_w));
               when 'n' => Add (Matrix, To_Symop (Translation_n));
               when '1' => 
                  case Rotation is 
                     when '3' => 
                        Add (Matrix, To_Symop (Translations_3_1, Axis_Number));
                     when '4' => 
                        Add (Matrix, To_Symop (Translations_4_1, Axis_Number));
                     when '6' => 
                        Add (Matrix, To_Symop (Translations_6_1, Axis_Number));
                     when others =>
                        raise UNKNOWN_ROTATION
                          with "mismatching translation " & 
                          Tr'Image & " for rotation " & Rotation'Image;
                  end case;
               when '2' => 
                  case Rotation is 
                     when '3' => 
                        Add (Matrix, To_Symop (Translations_3_2, Axis_Number));
                     when '6' => 
                        Add (Matrix, To_Symop (Translations_6_2, Axis_Number));
                     when others =>
                        raise UNKNOWN_ROTATION
                          with "mismatching translation " & 
                          Tr'Image & " for rotation " & Rotation'Image;
                  end case;
               when '3' => 
                  case Rotation is 
                     when '4' => 
                        Add (Matrix, To_Symop (Translations_4_3, Axis_Number));
                     when others =>
                        raise UNKNOWN_ROTATION
                          with "mismatching translation " & 
                          Tr'Image & " for rotation " & Rotation'Image;
                  end case;
               when '4' => 
                  case Rotation is 
                     when '6' => 
                        Add (Matrix, To_Symop (Translations_6_4, Axis_Number));
                     when others =>
                        raise UNKNOWN_ROTATION
                          with "mismatching translation " & 
                          Tr'Image & " for rotation " & Rotation'Image;
                  end case;
               when '5' => 
                  case Rotation is 
                     when '6' => 
                        Add (Matrix, To_Symop (Translations_6_5, Axis_Number));
                     when others =>
                        raise UNKNOWN_ROTATION
                          with "mismatching translation " & 
                          Tr'Image & " for rotation " & Rotation'Image;
                  end case;
               when others =>
                  raise UNKNOWN_TRANSLATION
                    with "translation character " & Tr'Image;
            end case;
         end loop;
         
      end;
      
      Inversion : Character;
      Rotation : Character;
      Axis : Character;
      Translations : String (1..2) := (others => ' ');
      
   begin
      Skip_Spaces (Symbol, Pos);
      Get_Inversion_Character (Inversion);
      Get_Rotation_Character (Rotation);
      Get_Axis_Character (Axis, Axis_Number);
      Get_Translation_Characters (Translations);
      
      if Rotation /= ' ' and then Axis /= ' ' then
         N_Rotations := N_Rotations + 1;
         Construct_Rotation_Matrix (Rotations (N_Rotations), 
                                    Inversion, Axis, Rotation,
                                    Translations, Preceeding_Axis);
      end if;
   end;
   
   function Decode_Hall (Symbol : in String) return Symop_Array is
      Max_Symops : constant Integer := 96;
      Symops : Symop_Array (1 .. Max_Symops);
      N_Symops : Positive := 1; -- The first element is allocated for the
                                -- unity matrix.
      Pos : Positive := 1;      -- current position in the string 'Symbol'.
      
      N_Inversions : Positive;
      
      Rotations : Symop_Array (1..3);
      N_Rotations : Natural := 0;
      
      Centering : Symop_Array (1..4);
      N_Centering : Positive;
      
      Preceeding_Axis : Natural := 0;
      
   begin
      Symops (1) := Unity_Matrix;
      
      Get_Hall_Symbol_Inversions (Symbol, Pos, N_Inversions);
      Get_Hall_Symbol_Centerings (Symbol, Pos, Centering, N_Centering);
      Get_Hall_Symbol_Rotations  (Symbol, Pos, Rotations, N_Rotations, Preceeding_Axis, 1);
      Get_Hall_Symbol_Rotations  (Symbol, Pos, Rotations, N_Rotations, Preceeding_Axis, 2);
      Get_Hall_Symbol_Rotations  (Symbol, Pos, Rotations, N_Rotations, Preceeding_Axis, 3);
      
      Put_Line ("Inversions:");
      for I in 1..N_Inversions loop
         Put (Inversion_Matrices (I));
         New_Line;
      end loop;
      
      Put_Line ("Centerings:");
      for I in 1..N_Centering loop
         Put (Centering (I));
         New_Line;
      end loop;
      
      Put_Line ("Rotations:");
      for I in 1..N_Rotations loop
         Put_Line (I'Image);
         Put (Rotations (I));
         New_Line;
      end loop;
      
      return Symops (1..N_Symops);
   end;
   
begin
   
   for I in 1 .. Argument_Count loop
      Put_Line (Argument (I));
      declare
         Symops : Symop_Array := Decode_Hall (Argument (I));
      begin
         Put_Line ("Symops:");
         for I in Symops'Range loop
            Put (Symops (I));
            New_Line;
         end loop;
      end;
   end loop;
   
end Decode_Hall;
