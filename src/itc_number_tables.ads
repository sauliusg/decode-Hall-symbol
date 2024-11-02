package ITC_Number_Tables is
   
   type ITC_Number_And_Cell_Symbol_Type is new String (1 .. 5);
   type Shmueli_Symbol_Type is new String (1 .. 24);
   
   type Cell_Symbol_Interpretation is record
      ITC_Number_And_Cell_Symbol : ITC_Number_And_Cell_Symbol_Type;
      Shmueli_Symbol : Shmueli_Symbol_Type;
   end record;
   
   function "+" (S : String) return ITC_Number_And_Cell_Symbol_Type is
      (
       ITC_Number_And_Cell_Symbol_Type
         (
          S (S'First .. S'Last) &
            (1 .. ITC_Number_And_Cell_Symbol_Type'Length - S'Length => ' ')
         )
      );
   
   function "+" (S : String) return Shmueli_Symbol_Type is
      (
       Shmueli_Symbol_Type
         (
          S (S'First .. S'Last) & 
            (1 .. Shmueli_Symbol_Type'Length - S'Length => ' ')
         )
      );
   
   type ITC_Number_And_Cell_Symbol_Array is array (Positive range <>)
     of Cell_Symbol_Interpretation;
   
   -- Default settings chosen for the 230 3D space groups; settings
   --  and cell choices are compatible with the International Tables
   --  of Crystallography, Vol. A:
   
   ITC_Number_Shmueli_Symbols : constant array (Positive range <>)
     of Shmueli_Symbol_Type :=
     (
      +"PAN$P1A000", -- space group No. 1
      +"PAC$I1A000", -- space group No. 2
      +"PMN$P2B000", -- space group No. 3
      +"PMN$P2B060", -- space group No. 4
      +"CMN$P2B000", -- space group No. 5
      +"PMN$I2B000", -- space group No. 6
      +"PMN$I2B006", -- space group No. 7
      +"CMN$I2B000", -- space group No. 8
      +"CMN$I2B006", -- space group No. 9
      +"PMC$I1A000$P2B000", -- space group No. 10
      +"PMC$I1A000$P2B060", -- space group No. 11
      +"CMC$I1A000$P2B000", -- space group No. 12
      +"PMC$I1A000$P2B006", -- space group No. 13
      +"PMC$I1A000$P2B066", -- space group No. 14
      +"CMC$I1A000$P2B006", -- space group No. 15
      +"PON$P2C000$P2A000", -- space group No. 16
      +"PON$P2C006$P2A000", -- space group No. 17
      +"PON$P2C000$P2A660", -- space group No. 18
      +"PON$P2C606$P2A660", -- space group No. 19
      +"CON$P2C006$P2A000", -- space group No. 20
      +"CON$P2C000$P2A000", -- space group No. 21
      +"FON$P2C000$P2A000", -- space group No. 22
      +"ION$P2C000$P2A000", -- space group No. 23
      +"ION$P2C606$P2A660", -- space group No. 24
      +"PON$P2C000$I2A000", -- space group No. 25
      +"PON$P2C006$I2A000", -- space group No. 26
      +"PON$P2C000$I2A006", -- space group No. 27
      +"PON$P2C000$I2A600", -- space group No. 28
      +"PON$P2C006$I2A606", -- space group No. 29
      +"PON$P2C000$I2A066", -- space group No. 30
      +"PON$P2C606$I2A000", -- space group No. 31
      +"PON$P2C000$I2A660", -- space group No. 32
      +"PON$P2C006$I2A666", -- space group No. 33
      +"PON$P2C000$I2A666", -- space group No. 34
      +"CON$P2C000$I2A000", -- space group No. 35
      +"CON$P2C006$I2A000", -- space group No. 36
      +"CON$P2C000$I2A006", -- space group No. 37
      +"AON$P2C000$I2A000", -- space group No. 38
      +"AON$P2C000$I2A060", -- space group No. 39
      +"AON$P2C000$I2A600", -- space group No. 40
      +"AON$P2C000$I2A660", -- space group No. 41
      +"FON$P2C000$I2A000", -- space group No. 42
      +"FON$P2C000$I2A333", -- space group No. 43
      +"ION$P2C000$I2A000", -- space group No. 44
      +"ION$P2C000$I2A660", -- space group No. 45
      +"ION$P2C000$I2A600", -- space group No. 46
      +"POC$I1A000$P2C000$P2A000", -- space group No. 47
      +"POC$I1A666$P2C000$P2A000", -- space group No. 48
      +"POC$I1A000$P2C000$P2A006", -- space group No. 49
      +"POC$I1A660$P2C000$P2A000", -- space group No. 50
      +"POC$I1A000$P2C600$P2A600", -- space group No. 51
      +"POC$I1A000$P2C600$P2A066", -- space group No. 52
      +"POC$I1A000$P2C606$P2A000", -- space group No. 53
      +"POC$I1A000$P2C600$P2A606", -- space group No. 54
      +"POC$I1A000$P2C000$P2A660", -- space group No. 55
      +"POC$I1A000$P2C660$P2A606", -- space group No. 56
      +"POC$I1A000$P2C006$P2A060", -- space group No. 57
      +"POC$I1A000$P2C000$P2A666", -- space group No. 58
      +"POC$I1A660$P2C000$P2A660", -- space group No. 59
      +"POC$I1A000$P2C666$P2A660", -- space group No. 60
      +"POC$I1A000$P2C606$P2A660", -- space group No. 61
      +"POC$I1A000$P2C606$P2A666", -- space group No. 62
      +"COC$I1A000$P2C006$P2A000", -- space group No. 63
      +"COC$I1A000$P2C066$P2A000", -- space group No. 64
      +"COC$I1A000$P2C000$P2A000", -- space group No. 65
      +"COC$I1A000$P2C000$P2A006", -- space group No. 66
      +"COC$I1A000$P2C060$P2A000", -- space group No. 67
      +"COC$I1A066$P2C660$P2A660", -- space group No. 68
      +"FOC$I1A000$P2C000$P2A000", -- space group No. 69
      +"FOC$I1A333$P2C000$P2A000", -- space group No. 70
      +"IOC$I1A000$P2C000$P2A000", -- space group No. 71
      +"IOC$I1A000$P2C000$P2A660", -- space group No. 72
      +"IOC$I1A000$P2C606$P2A660", -- space group No. 73
      +"IOC$I1A000$P2C060$P2A000", -- space group No. 74
      +"PTN$P4C000", -- space group No. 75
      +"PTN$P4C003", -- space group No. 76
      +"PTN$P4C006", -- space group No. 77
      +"PTN$P4C009", -- space group No. 78
      +"ITN$P4C000", -- space group No. 79
      +"ITN$P4C063", -- space group No. 80
      +"PTN$I4C000", -- space group No. 81
      +"ITN$I4C000", -- space group No. 82
      +"PTC$I1A000$P4C000", -- space group No. 83
      +"PTC$I1A000$P4C006", -- space group No. 84
      +"PTC$I1A660$P4C660", -- space group No. 85
      +"PTC$I1A666$P4C666", -- space group No. 86
      +"ITC$I1A000$P4C000", -- space group No. 87
      +"ITC$I1A063$P4C063", -- space group No. 88
      +"PTN$P4C000$P2A000", -- space group No. 89
      +"PTN$P4C660$P2A660", -- space group No. 90
      +"PTN$P4C003$P2A006", -- space group No. 91
      +"PTN$P4C663$P2A669", -- space group No. 92
      +"PTN$P4C006$P2A000", -- space group No. 93
      +"PTN$P4C666$P2A666", -- space group No. 94
      +"PTN$P4C009$P2A006", -- space group No. 95
      +"PTN$P4C669$P2A663", -- space group No. 96
      +"ITN$P4C000$P2A000", -- space group No. 97
      +"ITN$P4C063$P2A063", -- space group No. 98
      +"PTN$P4C000$I2A000", -- space group No. 99
      +"PTN$P4C000$I2A660", -- space group No. 100
      +"PTN$P4C006$I2A006", -- space group No. 101
      +"PTN$P4C666$I2A666", -- space group No. 102
      +"PTN$P4C000$I2A006", -- space group No. 103
      +"PTN$P4C000$I2A666", -- space group No. 104
      +"PTN$P4C006$I2A000", -- space group No. 105
      +"PTN$P4C006$I2A660", -- space group No. 106
      +"ITN$P4C000$I2A000", -- space group No. 107
      +"ITN$P4C000$I2A006", -- space group No. 108
      +"ITN$P4C063$I2A666", -- space group No. 109
      +"ITN$P4C063$I2A660", -- space group No. 110
      +"PTN$I4C000$P2A000", -- space group No. 111
      +"PTN$I4C000$P2A006", -- space group No. 112
      +"PTN$I4C000$P2A660", -- space group No. 113
      +"PTN$I4C000$P2A666", -- space group No. 114
      +"PTN$I4C000$P2D000", -- space group No. 115
      +"PTN$I4C000$P2D006", -- space group No. 116
      +"PTN$I4C000$P2D660", -- space group No. 117
      +"PTN$I4C000$P2D666", -- space group No. 118
      +"ITN$I4C000$P2D000", -- space group No. 119
      +"ITN$I4C000$P2D006", -- space group No. 120
      +"ITN$I4C000$P2A000", -- space group No. 121
      +"ITN$I4C000$P2A609", -- space group No. 122
      +"PTC$I1A000$P4C000$P2A000", -- space group No. 123
      +"PTC$I1A000$P4C000$P2A006", -- space group No. 124
      +"PTC$I1A660$P4C000$P2A000", -- space group No. 125
      +"PTC$I1A666$P4C000$P2A000", -- space group No. 126
      +"PTC$I1A000$P4C000$P2A660", -- space group No. 127
      +"PTC$I1A000$P4C000$P2A666", -- space group No. 128
      +"PTC$I1A660$P4C660$P2A660", -- space group No. 129
      +"PTC$I1A660$P4C660$P2A666", -- space group No. 130
      +"PTC$I1A000$P4C006$P2A000", -- space group No. 131
      +"PTC$I1A000$P4C006$P2A006", -- space group No. 132
      +"PTC$I1A666$P4C666$P2A006", -- space group No. 133
      +"PTC$I1A666$P4C666$P2A000", -- space group No. 134
      +"PTC$I1A000$P4C006$P2A660", -- space group No. 135
      +"PTC$I1A000$P4C666$P2A666", -- space group No. 136
      +"PTC$I1A666$P4C666$P2A666", -- space group No. 137
      +"PTC$I1A666$P4C666$P2A660", -- space group No. 138
      +"ITC$I1A000$P4C000$P2A000", -- space group No. 139
      +"ITC$I1A000$P4C000$P2A006", -- space group No. 140
      +"ITC$I1A063$P4C063$P2A063", -- space group No. 141
      +"ITC$I1A063$P4C063$P2A069", -- space group No. 142
      +"PRN$P3C000", -- space group No. 143
      +"PRN$P3C004", -- space group No. 144
      +"PRN$P3C008", -- space group No. 145
      +"RRN$P3C000", -- space group No. 146
      +"PRC$I3C000", -- space group No. 147
      +"RRC$I3C000", -- space group No. 148
      +"PRN$P3C000$P2G000", -- space group No. 149
      +"PRN$P3C000$P2F000", -- space group No. 150
      +"PRN$P3C004$P2G000", -- space group No. 151
      +"PRN$P3C004$P2F008", -- space group No. 152
      +"PRN$P3C008$P2G000", -- space group No. 153
      +"PRN$P3C008$P2F004", -- space group No. 154
      +"RRN$P3C000$P2F000", -- space group No. 155
      +"PRN$P3C000$I2F000", -- space group No. 156
      +"PRN$P3C000$I2G000", -- space group No. 157
      +"PRN$P3C000$I2F006", -- space group No. 158
      +"PRN$P3C000$I2G006", -- space group No. 159
      +"RRN$P3C000$I2F000", -- space group No. 160
      +"RRN$P3C000$I2F006", -- space group No. 161
      +"PRC$I3C000$P2G000", -- space group No. 162
      +"PRC$I3C000$P2G006", -- space group No. 163
      +"PRC$I3C000$P2F000", -- space group No. 164
      +"PRC$I3C000$P2F006", -- space group No. 165
      +"RRC$I3C000$P2F000", -- space group No. 166
      +"RRC$I3C000$P2F006", -- space group No. 167
      +"PHN$P6C000", -- space group No. 168
      +"PHN$P6C002", -- space group No. 169
      +"PHN$P6C005", -- space group No. 170
      +"PHN$P6C004", -- space group No. 171
      +"PHN$P6C008", -- space group No. 172
      +"PHN$P6C006", -- space group No. 173
      +"PHN$I6C000", -- space group No. 174
      +"PHC$I1A000$P6C000", -- space group No. 175
      +"PHC$I1A000$P6C006", -- space group No. 176
      +"PHN$P6C000$P2F000", -- space group No. 177
      +"PHN$P6C002$P2F000", -- space group No. 178
      +"PHN$P6C005$P2F000", -- space group No. 179
      +"PHN$P6C004$P2F000", -- space group No. 180
      +"PHN$P6C008$P2F000", -- space group No. 181
      +"PHN$P6C006$P2F000", -- space group No. 182
      +"PHN$P6C000$I2F000", -- space group No. 183
      +"PHN$P6C000$I2F006", -- space group No. 184
      +"PHN$P6C006$I2F006", -- space group No. 185
      +"PHN$P6C006$I2F000", -- space group No. 186
      +"PHN$I6C000$P2G000", -- space group No. 187
      +"PHN$I6C006$P2G000", -- space group No. 188
      +"PHN$I6C000$P2F000", -- space group No. 189
      +"PHN$I6C006$P2F000", -- space group No. 190
      +"PHC$I1A000$P6C000$P2F000", -- space group No. 191
      +"PHC$I1A000$P6C000$P2F006", -- space group No. 192
      +"PHC$I1A000$P6C006$P2F006", -- space group No. 193
      +"PHC$I1A000$P6C006$P2F000", -- space group No. 194
      +"PCN$P3Q000$P2C000$P2A000", -- space group No. 195
      +"FCN$P3Q000$P2C000$P2A000", -- space group No. 196
      +"ICN$P3Q000$P2C000$P2A000", -- space group No. 197
      +"PCN$P3Q000$P2C606$P2A660", -- space group No. 198
      +"ICN$P3Q000$P2C606$P2A660", -- space group No. 199
      +"PCC$I3Q000$P2C000$P2A000", -- space group No. 200
      +"PCC$I3Q666$P2C000$P2A000", -- space group No. 201
      +"FCC$I3Q000$P2C000$P2A000", -- space group No. 202
      +"FCC$I3Q333$P2C000$P2A000", -- space group No. 203
      +"ICC$I3Q000$P2C000$P2A000", -- space group No. 204
      +"PCC$I3Q000$P2C606$P2A660", -- space group No. 205
      +"ICC$I3Q000$P2C606$P2A660", -- space group No. 206
      +"PCN$P3Q000$P4C000$P2D000", -- space group No. 207
      +"PCN$P3Q000$P4C666$P2D666", -- space group No. 208
      +"FCN$P3Q000$P4C000$P2D000", -- space group No. 209
      +"FCN$P3Q000$P4C993$P2D939", -- space group No. 210
      +"ICN$P3Q000$P4C000$P2D000", -- space group No. 211
      +"PCN$P3Q000$P4C939$P2D399", -- space group No. 212
      +"PCN$P3Q000$P4C393$P2D933", -- space group No. 213
      +"ICN$P3Q000$P4C393$P2D933", -- space group No. 214
      +"PCN$P3Q000$I4C000$I2D000", -- space group No. 215
      +"FCN$P3Q000$I4C000$I2D000", -- space group No. 216
      +"ICN$P3Q000$I4C000$I2D000", -- space group No. 217
      +"PCN$P3Q000$I4C666$I2D666", -- space group No. 218
      +"FCN$P3Q000$I4C666$I2D666", -- space group No. 219
      +"ICN$P3Q000$I4C939$I2D399", -- space group No. 220
      +"PCC$I3Q000$P4C000$P2D000", -- space group No. 221
      +"PCC$I3Q666$P4C000$P2D000", -- space group No. 222
      +"PCC$I3Q000$P4C666$P2D666", -- space group No. 223
      +"PCC$I3Q666$P4C666$P2D666", -- space group No. 224
      +"FCC$I3Q000$P4C000$P2D000", -- space group No. 225
      +"FCC$I3Q000$P4C666$P2D666", -- space group No. 226
      +"FCC$I3Q333$P4C993$P2D939", -- space group No. 227
      +"FCC$I3Q999$P4C993$P2D939", -- space group No. 228
      +"ICC$I3Q000$P4C000$P2D000", -- space group No. 229
      +"ICC$I3Q000$P4C393$P2D933"  -- space group No. 230
     );
   
   -- Alternative cell coices for some 3D space groups:
   
   ITC_Number_And_Cell_Symbol_Table : constant 
     ITC_Number_And_Cell_Symbol_Array :=
     (
      (
       +"1:2",
       +"BAN$P1A000"
      ),
      (
       +"1:3",
       +"CAN$P1A000"
      ),
      (
       +"1:4",
       +"FAN$P1A000"
      ),
      (
       +"1:5",
       +"IAN$P1A000"
      ),
      (
       +"2:2",
       +"BAN$I1A000"
      ),
      (
       +"2:3",
       +"CAN$I1A000"
      ),
      (
       +"2:4",
       +"FAN$I1A000"
      ),
      (
       +"2:5",
       +"IAN$I1A000"
      ),
      (
       +"48:2",
       +"POC$I1A000$P2C660$P2A066"
      ),
      (
       +"50:2",
       +"POC$I1A000$P2C660$P2A060"
      ),
      (
       +"59:2",
       +"POC$I1A000$P2C660$P2A600"
      ),
      (
       +"68:2",
       +"COC$I1A000$P2C600$P2A606"
      ),
      (
       +"70:2",
       +"FOC$I1A000$P2C990$P2A099"
      ),
      (
       +"85:2",
       +"PTC$I1A000$P4C600"
      ),
      (
       +"86:2",
       +"PTC$I1A000$P4C066"
      ),
      (
       +"88:2",
       +"ITC$I1A000$P4C933"
      ),
      (
       +"125:2",
       +"PTC$I1A000$P4C600$P2A060"
      ),
      (
       +"126:2",
       +"PTC$I1A000$P4C600$P2A066"
      ),
      (
       +"129:2",
       +"PTC$I1A000$P4C600$P2A600"
      ),
      (
       +"130:2",
       +"PTC$I1A000$P4C600$P2A606"
      ),
      (
       +"133:2",
       +"PTC$I1A000$P4C606$P2A060"
      ),
      (
       +"134:2",
       +"PTC$I1A000$P4C606$P2A066"
      ),
      (
       +"137:2",
       +"PTC$I1A000$P4C606$P2A600"
      ),
      (
       +"138:2",
       +"PTC$I1A000$P4C606$P2A606"
      ),
      (
       +"141:2",
       +"ITC$I1A000$P4C393$P2A000"
      ),
      (
       +"142:2",
       +"ITC$I1A000$P4C393$P2A006"
      ),
      (
       +"146:H",
       +"RRN$P3C000"
      ),
      (
       +"146:R",
       +"PRN$P3Q000"
      ),
      (
       +"148:H",
       +"RRC$I3C000"
      ),
      (
       +"148:R",
       +"PRC$I3Q000"
      ),
      (
       +"155:H",
       +"RRN$P3C000$P2F000"
      ),
      (
       +"155:R",
       +"PRN$P3Q000$P2E000"
      ),
      (
       +"160:H",
       +"RRN$P3C000$I2F000"
      ),
      (
       +"160:R",
       +"PRN$P3Q000$I2E000"
      ),
      (
       +"161:H",
       +"RRN$P3C000$I2F006"
      ),
      (
       +"161:R",
       +"PRN$P3Q000$I2E666"
      ),
      (
       +"166:H",
       +"RRC$I3C000$P2F000"
      ),
      (
       +"166:R",
       +"PRC$I3Q000$P2E000"
      ),
      (
       +"167:H",
       +"RRC$I3C000$P2F006"
      ),
      (
       +"167:R",
       +"PRC$I3Q000$P2E666"
      ),
      (
       +"201:2",
       +"PCC$I3Q000$P2C660$P2A066"
      ),
      (
       +"203:2",
       +"FCC$I3Q000$P2C330$P2A033"
      ),
      (
       +"222:2",
       +"PCC$I3Q000$P4C600$P2D006"
      ),
      (
       +"224:2",
       +"PCC$I3Q000$P4C066$P2D660"
      ),
      (
       +"227:2",
       +"FCC$I3Q000$P4C693$P2D936"
      ),
      (
       +"228:2",
       +"FCC$I3Q000$P4C093$P2D930"
      )
     );
   
end ITC_Number_Tables;
