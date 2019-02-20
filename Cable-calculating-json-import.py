#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      lavrinenko.i
#
# Created:     28.08.2018
# Copyright:   (c) lavrinenko.i 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import math
from scipy.special import cbrt
import os.path
import xlrd, xlwt, openpyxl
from math import radians, cos, sin, asin, sqrt
from PIL import Image, ImageDraw, ImageFont
from collections import Counter
import json
import pathlib

font = ImageFont.truetype("arial.ttf", 20)

PATH_FILE_JSON_1 = 'STB.json'
#POSELOK_DIR_PATH = PATH_FILE_JSON_1[:-5]
POSELOK_DIR_PATH = 'stolbishi'
##DIR_PATH = ''
#ANKERNYE_FILE = 'ankernye.txt'

IMAGE_FILE_PODDER = 'template-podderjivaushaya.jpg'
IMAGE_FILE_ANKER = 'template-ankernaya.jpg'

#IMAGE_SAVE_PATH = 'opori/'

IMAGE_SAVE_PATH = 'poselki/'+ POSELOK_DIR_PATH +'/calculated-poles/'

pathlib.Path(IMAGE_SAVE_PATH).mkdir(parents=True, exist_ok=True)


op_massive = []
trass_massive = []
opor_massive = []
leng_massive = []
ank_massive = []

#-------------------------------------------------------------------------------

g = 9.8 # м/с^2
m = 68.6 # кг/км  - Масса кабеля, кг/км
Sp = 2      # стрела провеса в % от длины пролета. Стрела провеса начальная, % от длины пролета
h = 0   #   Перепад высот при максимальной длине пролета, м
E_kab = 4.76 * 10**3 #   Модуль упругости начальный, Н/мм2
D_kab = 9.3  # Наружный диаметр кабеля, мм
TKLR_d = 19.82 * 10**-6 # Температурный коэффициент линейного расширения, *10 ^-6, 1/°С
H_z = 5.5 # Высота подвеса кабеля

#-------------------------------------------------------------------------------

def Sort_by_Op(stri):
    stri = stri[0][6:]
    return int(stri)

def toFixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"

def TKLR(T, T_sr):
        tklr = (T - T_sr)*TKLR_d
        #print('2.6.1 Температурный коэффициент расширения:',tklr)
        return(tklr)

def haversine(lat1, lon1, lat2, lon2):
    """
    Вычисляет расстояние в километрах между двумя точками, учитывая окружность Земли.
    https://en.wikipedia.org/wiki/Haversine_formula
    """

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, (lon1, lat1, lon2, lat2))

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    km = 6367 * c * 1000
    return km

##def create_massive_from_xls(PATH_FILE_XLS):
##    i = 0
##    massive = []
##    rb = xlrd.open_workbook(PATH_FILE_XLS)
##    #print(rb)
##    wb = openpyxl.load_workbook(filename = PATH_FILE_XLS)
##    sheet = wb.active
##    i = 3
##    for row in sheet.rows:
##        while sheet['A'+str(i)].value != None:
##            lat = str(sheet['A'+str(i)].value)
##            lon = str(sheet['B'+str(i)].value)
##            holder_id = sheet['E'+str(i)].value
##            massive.append([bin(holder_id)[2:], lat, lon])
##            i = i + 1
##    return massive

##def create_massive_of_lines(PATH_FILE_JSON):
##    mas = []
##    for line in open(PATH_FILE_JSON):
##        mas.append(line)
##        print(line)
##        first = line.split(',')[1][5:-1]
##        print(first)
##        second = line.split(',')[3][5:-2]
##        print(second)
##        #subline = [first, second]
##        trass_massive.append([ 'trassa-' + str(int(line.split(',')[1][5:-1],2)) + '-' + str(int(line.split(',')[3][5:-2],2)) , float(toFixed(haversine(float(find_locate(first)[0]), float(find_locate(first)[1]), float(find_locate(second)[0]), float(find_locate(second)[1]) ), 2)), line.split(',')[1][5:-1], line.split(',')[3][5:-2] ])
##        #print(trass_massive)
##    return trass_massive

def create_massive_of_lines_1(PATH_FILE_JSON):
    mas = []
    for line in open(PATH_FILE_JSON):
        data = json.loads(line)
        for obj in data["objects"]:
            if obj['type'] == 301:
                lat = obj['location']['latitude']
                lon = obj['location']['longitude']
                mas.append( [ obj['id'],
                              lat,
                              lon,
                              obj['number'] ])
                #print(obj['id'], lat, lon)
    #print(mas)
    return mas

def id_to_nom(ide):
    #print(ide)
    #print(op_massive)
    for it in op_massive:
        #print(it)
        if ide == it[0]:
            return(it[3])


def create_massive_of_lines_2(PATH_FILE_JSON):
    mas = []
    for line in open(PATH_FILE_JSON):
        data = json.loads(line)
        for obj in data["objects"]:
            if obj['type'] == 401:
                #print(obj['first']['id'],obj['second']['id'])
                print(obj['first']['id'])
                first = id_to_nom(obj['first']['id'])
                #print(first)
                #print(op_massive)
                print(first)
                second = id_to_nom(obj['second']['id'])
                print(second)
                #subline = [first, second]
                trass_massive.append([
                                        'trassa-' + str(first) + '-' + str(second) ,
                                         float(toFixed(haversine(float(find_locate(str(first))[0]) ,
                                                                 float(find_locate(first)[1]),
                                                                 float(find_locate(second)[0]),
                                                                 float(find_locate(second)[1]))
                                                , 2)),
                                        obj['first']['id'],
                                        obj['second']['id']
                                    ])
                print(trass_massive)
    return trass_massive


def create_massive_of_ankers(PATH_FILE):
    mas_ankers = []
    for line in open(PATH_FILE):
        mas_ankers.append(line[:-1])
        print(line)
    return mas_ankers

##def create_massive_of_ends(trass_massive):


def L_ot_T(T_montajnaja):
    L_montajnaja = L_n0 * (1 + TKLR(T_montajnaja,T_sr))
    return(L_montajnaja)

def Provis(L_mot):
    a = 3 * (L**2 + (h**2 / 2) - (L * L_mot))  / 8
    b = (-3 * W * L**3 * L_mot) / (64 * E_kab * S_kab)
    ed = (a/3)**3 + (-b/2)**2
    if ed < 0:
        S_mo = 2 * math.sqrt((-1 * a)/3) * math.cos( (1/3) / math.cos((-1 * b/2)/(-1 * a/3)**1.500))
    else:
        sqr = (a/3)**3 + (-1 * b/2)**2
        S_mo = cbrt((-1 * b/2) + math.sqrt(sqr)) + cbrt((-1 * b/2) - math.sqrt(sqr))
    return(S_mo)

def Nagruzka(S_mo):
        H_i = W * (L**2) / (8 * S_mo)
        return(H_i)

def Nagruzka_ot_dlini(S_mo, L):
        H_j = W * (L**2) / (8 * S_mo)
        return(H_j)

def find_locate(id_s):
    for item in op_massive:
        if str(id_s) == str(item[3]):
            return item[1] , item[2]

def find_lingth(stolb_number):
    leng_mas = trass_massive
    trass_st = []
    for trass in leng_mas:
        if (id_to_nom(trass[2]) == stolb_number) or (id_to_nom(trass[3]) == stolb_number):
            trass_st.append(trass)
    return trass_st

def opor_massive_create(trass_massive):
    print(trass_massive)
    cut_op = []
    for op in trass_massive:
        cut_op.append(id_to_nom(op[2]))
        cut_op.append(id_to_nom(op[3]))
    print('ccc', cut_op)
    c = Counter(cut_op)
    print(list(c))
    print(c)
    print(id_to_nom(trass_massive[0][2]))
    op_mas_nag = []
    re = 0
    for sto in list(c):
        print(sto)
        print(c[sto])

        if c[sto] == 1:
            t = find_lingth(sto)
            for one in t:
                op_mas_nag.append(['opora-'+ str(sto) , one[1] / 2 , one[1] , 0 ])
                if str(sto) not in ank_massive:
                    ank_massive.append(str(sto))
                    print('dob',ank_massive)
        elif c[sto] == 2:
            t = find_lingth(sto)
            op_mas_nag.append(['opora-'+ str(sto) , (t[0][1] / 2 + t[1][1] /2 ), t[0][1] , t[1][1]])
        elif c[sto] == 3:
            t = find_lingth(sto)
            op_mas_nag.append(['opora-'+ str(sto) , (t[0][1] / 2 + t[1][1] /2 ), t[0][1] , t[1][1]])
            op_mas_nag.append(['opora-'+ str(sto) , t[2][1] / 2 , t[2][1] , 0 ])
            if str(sto) not in ank_massive:
                    ank_massive.append(str(sto))
        elif c[sto] == 4:
            t = find_lingth(sto)
            op_mas_nag.append(['opora-'+ str(sto) , (t[0][1] / 2 + t[1][1] /2 ), t[0][1] , t[1][1]])
            op_mas_nag.append(['opora-'+ str(sto) , (t[2][1] / 2 + t[3][1] /2 ) , t[2][1] , t[3][1] ])
            if str(sto) not in ank_massive:
                    ank_massive.append(str(sto))
            #op_mas_nag.append(['opora-'+ str(sto) , float(toFixed(trass_massive[(re)][4] / 2  + trass_massive[(re + 1)][4] / 2 , 2 )) , float(trass_massive[(re)][4]) , float(trass_massive[(re + 1)][4]) ])
##    op_mas_nag = []
##    re = 0
##    while (re + 1) < len(trass_massive):
##        if trass_massive[re][3] == trass_massive[(re + 1)][2]:
##            op_mas_nag.append(['opora-'+ str(int(trass_massive[(re + 1)][2] , 2)) , float(toFixed(trass_massive[(re)][4] / 2  + trass_massive[(re + 1)][4] / 2 , 2 )) , float(trass_massive[(re)][4]) , float(trass_massive[(re + 1)][4]) ])
##        else:
##            op_mas_nag.append(['opora-'+ str(int(trass_massive[(re)][2] , 2)) , float(toFixed(trass_massive[(re)][1] / 2 , 2 )) , float(trass_massive[(re)][1]) , 0 ])
##        re = re + 1
    return op_mas_nag

def moment_sili(F):
    mom_mas = []
    for pe in F:
        mom_mas.append( toFixed((float(pe) * H_z) , 2) )
    return mom_mas


#create_massive_of_lines_1(PATH_FILE_JSON_1)

#op_massive = create_massive_from_xls(PATH_FILE_XLS)

op_massive = create_massive_of_lines_1(PATH_FILE_JSON_1)

#ank_massive = create_massive_of_ankers(DIR_PATH + ANKERNYE_FILE)
print(ank_massive)
##print('first', op_massive)
create_massive_of_lines_2(PATH_FILE_JSON_1)

print('Трассы Кузиметьево')
for j in trass_massive:
    print(j[0] , j[1])

## 2.1 Вес кабеля - нагрузка каждого погонного метра  W

W = ( m * g )/1000
print('2.1 Вес кабеля:', W,'Н/м')

for it in trass_massive:
    L = it[1]
    print('----------------Рассчет трассы '+str(it[0])+' длинной '+str(it[1])+'м ------------------')

    ## 2.2. Растягивающая нагрузка, действующая на кабель H
    S = (L * Sp/100)
    H = W * (L**2) / (8 * S)
    print ('2.2. Растягивающая нагрузка, действующая на кабель:', H ,'Н')

    ## 2.3. Перепад высот между опорами
    H_first = H
    L1 = L - ( 2 * h * H_first / (W * L))
    L2 = L - ( 2 * h * H_first / (W * L))
    print('Длина нагружающего пролета с каждой из сторон', L1/2)
    S1 = W * L1**2 / (8 * H)
    S2 = W * L2**2 / (8 * H)
    print('2.3 Перепад высот между опорами:',S1, S2)

    ## 2.4. Длина подвешенного кабеля.
    L_kab = L + (4/3 * ((S1**2/L1)+(S2**2/L2)))
    print('2.4. Длина подвешенного кабеля:',L_kab, 'м')

    ## 2.5. Длина кабеля в ненагруженном состоянии:

    S_kab = math.pi * (D_kab/2)**2  # вычисляем площадь кабеля в мм2
    L_n0 = L_kab / (1 + (H/(E_kab * S_kab)))
    print('2.5. Длина кабеля в ненагруженном состоянии:', L_n0, 'м')

    ## 2.6. Длина кабеля в ненагруженном состоянии с учетом температуры:
    T = 10
    T_sr = 10  # Среднегодовая температура в Казани  4.6
    L_nh = L_n0 * (1 + TKLR(T,T_sr))
    print('2.6. Длина кабеля в ненагруженном состоянии с учетом температуры:', L_nh, 'м')
    it.append(L_nh)


    ## 2.7. Вес кабеля при воздействии максимального гололеда
    Ro_l = 900 # плотность льда
    Ci = 15  # толщина стенки льда в мм
    K_i = 1
    K_d = 1
    G_gol = Ro_l * g * math.pi * K_i * K_d * Ci * (D_kab + K_i * K_d * Ci)/10**6    # Вес льда
    W_gol = W + Ro_l * g * math.pi * K_i * K_d * Ci * (D_kab + K_i * K_d * Ci)/10**6
    print('2.7. Вес кабеля при воздействии максимального гололеда:',W_gol, 'Н')

    ##  2.8. Ветровая нагрузка на кабель при гололеде
    a_w = 0.71
    K_l = 1.2  # коэффициент, учитывающий влияние длины пролета на ветровую нагрузку, равный 1,2 при длине пролета до 50 м, 1,1 - при 100 м, 1,05 - при 150 м, 1,0 - при 250 м и более (промежуточные значения Kl определяются интерполяцией);
    K_w = 1
    C_x = 1.2
    V_vetra = 29
    angles = [ 45 , 90 ]
    W_vetr_mas = []



    W_vetr_davlenie1 = V_vetra**2 / 1.6
    W_vetr_davlenie = 500


    for angle in angles:

        W_vetr = a_w * K_l * K_w * C_x * W_vetr_davlenie * D_kab * (math.sin(math.pi * angle / 180 ))**2 / 10**3
        W_vetr_mas.append(W_vetr)




    W_vetr_ygol = a_w * K_l * K_w * C_x * 200 * (D_kab + 2 * K_i * K_d * Ci) * (math.sin(math.pi * 90 / 180 ))**2 / 10**3


    print('2.8. Ветровая нагрузка на кабель при гололеде:', W_vetr_ygol, 'Н/м')


    #print('2.8. Ветровая нагрузка на кабель при угле', str(angle)+'°:', W_vetr_ygol, 'Н/м')

    ##  2.9. Максимальная нагрузка, действующая на кабель.

    W_max = math.sqrt(W_gol**2 + W_vetr_ygol**2)
    print('2.9. Максимальная нагрузка, действующая на кабель:', W_max, 'Н/м')


    ##  2.10 Расчет максимальной стрелы провеса.

    ## 2.10.1 Определив максимальную нагрузку, можно узнать длину кабеля в нагруженном состоянии (из п. 2.5. и п. 2.2.):

    a = 3 * (L**2 + (h**2 / 2) - (L * L_nh))  / 8
    b = (-3 * W_max * L**3 * L_nh) / (64 * E_kab * S_kab)

    ed = (a/3)**3 + (-b/2)**2


    if ed < 0:
        S_max = 2 * math.sqrt(-a/3) * math.cos( (1/3) / math.cos((-b/2)/(-a/3)**1.5))
        print('2.10 Расчет максимальной стрелы провеса:',S_max, 'м')
    else:
        S_max = ((-b/2) + math.sqrt((a/3)**3 + (-b/2)**2))**(1/3.0) + ((-b/2) - math.sqrt((a/3)**3 + (-b/2)**2))**(1/3.0)
        print('2.10 Расчет максимальной стрелы провеса:',S_max, 'м')


    ##  2.11 Максимальная растягивающая нагрузка при наихудших условиях.

    H_max = (W_max * L**2)/(8 * S_max)
    print('2.11 Максимальная растягивающая нагрузка при наихудших условия:',H_max, 'Н')


    ##  2.12. Расчет монтажной стрелы провеса, нагрузки и монтажной таблицы

    temp_list = [-60, -30, -20, -10, 0, 10, 20, 30, 40,	50,	60,	70]
    print("провис в гололед:", S_max,  )
    print()
    print('T , Расстояние от земли, Длина кабеля, Растягивающая нагрузка кабеля')
    for i in temp_list:
        print(str(i)+'°C:', toFixed(( H_z - float(Provis(L_ot_T(i)))), 2) , toFixed(L_ot_T(i), 2), toFixed(Nagruzka_ot_dlini(Provis(L_ot_T(i)), L), 2) )

##    M_vertic_nagr = H_max * H_z
##
##    print('Момент на основание опоры:', toFixed(M_vertic_nagr, 2), 'Н*м')
    print()


leng_massive = opor_massive_create(trass_massive)

Y_nw =1 # коэффициент надежности по ответственности линии равный 1 для линий до 220 кВ;
Y_p = 1 # региональный коэффициент, равный 1.
Y_f_v = [1.3,1.1]
Y_f_n = [1.3,1]
Y_f_g = 1.3
Y_d = [1,0.5]

print('fifi',trass_massive)
print(op_massive)
print(leng_massive)
print()

leng_massive.sort(key=Sort_by_Op)
print(leng_massive)
file_count = 1
for zapis in leng_massive:
    P_vetr_45_I_L = 0
    P_vetr_45_II_L = 0
    P_vetr_45_I_R = 0
    P_vetr_45_II_R = 0
    P_vetr_90_I_L = 0
    P_vetr_90_II_L = 0
    P_vetr_90_I_R = 0
    P_vetr_90_II_R = 0


    if str(zapis[0][6:]) not in ank_massive:
        image = Image.open(IMAGE_FILE_PODDER)
        print('----------------------------Расчет для опоры '+zapis[0]+'----------------------------')
        print('Длина сегмента слева '+toFixed(zapis[2] , 2) +' м, справа '+ toFixed(zapis[3] , 2  )+' м')
        print(zapis[0][6:])

        draw = ImageDraw.Draw(image)
        draw.text((1013, 177),zapis[0][6:],(0,0,0),font=font)                                                    # Печатает номер опоры
        draw.text((228, 408),toFixed(zapis[2] , 2),(0,0,0),font=font)                                            # Печатает длину пролета слева
        draw.text((592, 408),toFixed(zapis[3] , 2  ),(0,0,0),font=font)                                          # Печатает длину пролета справа



        draw.text((380, 525),'ДОТс-П-3кН',(0,0,0),font=font)                                                     # Печатает марку кабеля



        ## Рассчет для групп предельных состояний

        ## Идрисовская поебота 1

        G_n = m * g  * zapis[1] / 1000

        G_nI_II = []
        for Y_d_i in Y_d:
            G_nI_II.append( toFixed((G_n * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1

        #print(toFixed(G_n,2))
        print('Нормальный режим : Горизонтальная поперечная нагрузка, Н', [ 0, 0] )
        draw.text((319, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((701, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Нормальный режим : Вертикальная нагрузка, Н             ', G_nI_II )
        draw.text((412, 950),G_nI_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((806, 950),G_nI_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((412, 1050),G_nI_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((806, 1050),G_nI_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((412, 1145),G_nI_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((806, 1145),G_nI_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((412, 1312),G_nI_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((806, 1312),G_nI_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Нормальный режим : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
        draw.text((513, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((925, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((513, 1050),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((925, 1050),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((513, 1145),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((925, 1145),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((513, 1220),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((925, 1220),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((513, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((925, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((319, 1312),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((319, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((701, 1312),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((701, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((608, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((1026, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((608, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((1026, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Нормальный режим : Момент силы на основание опоры, Н    ', [ 0, 0] )
        print('----------------------------------------------------------------------------------')



        fd = 0

        for yg in W_vetr_mas:
            P_vetr_ygol = yg * zapis[1]

            P_vetr_ygol_II = []
            for Y_d_i in Y_d:
                P_vetr_ygol_II.append( toFixed((P_vetr_ygol * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1

            #print(toFixed(P_vetr_ygol,2))
            mom_mas = []
            if fd == 0:
                print('Режим максимального ветра под углом 45 градусов к линии : Горизонтальная поперечная нагрузка, Н' , P_vetr_ygol_II)
                draw.text((319, 1050),P_vetr_ygol_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((701, 1050),P_vetr_ygol_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('Режим максимального ветра под углом 45 градусов к линии : Вертикальная нагрузка, Н             ', G_nI_II )
                print('Режим максимального ветра под углом 45 градусов к линии : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
                print('Режим максимального ветра под углом 45 градусов к линии : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_II) )
                draw.text((608, 1050),moment_sili(P_vetr_ygol_II)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((1026, 1050),moment_sili(P_vetr_ygol_II)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('----------------------------------------------------------------------------------')
            elif fd == 1:
                print('Режим максимального ветра перпендикулярного линии : Горизонтальная поперечная нагрузка, Н' , P_vetr_ygol_II)
                draw.text((319, 1145),P_vetr_ygol_II[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((701, 1145),P_vetr_ygol_II[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('Режим максимального ветра перпендикулярного линии : Вертикальная нагрузка, Н             ', G_nI_II )
                print('Режим максимального ветра перпендикулярного линии : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
                print('Режим максимального ветра перпендикулярного линии : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_II) )
                draw.text((608, 1145),moment_sili(P_vetr_ygol_II)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((1026, 1145),moment_sili(P_vetr_ygol_II)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('----------------------------------------------------------------------------------')
            fd = fd + 1

        ## Расчет нагрузкина опору режим гололеда с ветром
        ## вертикальная нагрузка:

        W_vertikalnaya_gololed =  W_gol * zapis[1]
        W_vertikalnaya_gololed_mas = []
        for Y_d_i in Y_d:
            W_vertikalnaya_gololed_mas.append( toFixed((W_vertikalnaya_gololed * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1



        ## Горизонтальная нагрузка в гололед:

        P_vetr_ygol_gol = W_vetr_ygol * zapis[1]
        P_vetr_ygol_gol_mas = []
        for Y_d_i in Y_d:
            P_vetr_ygol_gol_mas.append( toFixed((P_vetr_ygol_gol * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1

        print('Режим гололеда с ветром : Горизонтальная поперечная нагрузка, Н', P_vetr_ygol_gol_mas)
        draw.text((319, 1220),P_vetr_ygol_gol_mas[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((701, 1220),P_vetr_ygol_gol_mas[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Режим гололеда с ветром : Вертикальная нагрузка, Н             ', W_vertikalnaya_gololed_mas )
        draw.text((412, 1220),W_vertikalnaya_gololed_mas[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((806, 1220),W_vertikalnaya_gololed_mas[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Режим гололеда с ветром : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
        print('Режим гололеда с ветром : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_gol_mas) )
        draw.text((608, 1220),moment_sili(P_vetr_ygol_gol_mas)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((1026, 1220),moment_sili(P_vetr_ygol_gol_mas)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('----------------------------------------------------------------------------------')


        ## Горизонтальная продольная нагрузка при обрыве    ***** нужно вычислять исходя из максимального пролета для данного столба



        d = []

        if zapis[2] >= zapis[3]:
            #print('Максимальная длина пролета для этого столба', zapis[2])
            inc = 0
            for Y_d_i in Y_d:
                d.append( toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[2] ) *  Y_d_i * Y_f_v[inc]), 2)  )
                draw.text((402, 467),toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[2] )), 2),(0,0,0),font=font)                                          # Печатает тяжение
                inc = inc + 1
        else:
            inc = 0
            for Y_d_i in Y_d:
                d.append( toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[3] ) *  Y_d_i * Y_f_v[inc]), 2)  )
                draw.text((402, 467),toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[3] )), 2),(0,0,0),font=font)                                          # Печатает тяжение
                inc = inc + 1
            #print('Максимальная длина пролета для этого столба', zapis[3])

        print('Режим одноосного тяжения (обрыв) : Горизонтальная поперечная нагрузка, Н', [ 0, 0])
        print('Режим одноосного тяжения (обрыв) : Вертикальная нагрузка, Н             ', G_nI_II )
        print('Режим одноосного тяжения (обрыв) : Горизонтальная продольная нагрузка, Н', d )
        draw.text((513, 1312),d[0],(0,0,0),font=font)
        draw.text((925, 1312),d[1],(0,0,0),font=font)
        print('Режим одноосного тяжения (обрыв) : Момент силы на основание опоры, Н    ', moment_sili(d) )
        draw.text((608, 1312),moment_sili(d)[0],(0,0,0),font=font)
        draw.text((1026, 1312),moment_sili(d)[1],(0,0,0),font=font)
        print('----------------------------------------------------------------------------------')



        ## Вес монтажника с оснасткой


        M_montajnika = 100
        P_montajnika = M_montajnika * 9.8

        print('Монтажный режим (вес монтажника с оснасткой) : Горизонтальная поперечная нагрузка, Н', [ 0, 0])
        print('Монтажный режим (вес монтажника с оснасткой) : Вертикальная нагрузка, Н             ', [ toFixed(P_montajnika + 24.24 , 2), toFixed(P_montajnika + 12.12, 2) ] )
        draw.text((402, 1380),toFixed(P_montajnika + float(G_nI_II[0]) , 2),(0,0,0),font=font)
        draw.text((806, 1380),toFixed(P_montajnika + float(G_nI_II[1]) , 2),(0,0,0),font=font)
        print('Монтажный режим (вес монтажника с оснасткой) : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
        print('Монтажный режим (вес монтажника с оснасткой) : Момент силы на основание опоры, Н    ', [ 0, 0] )
        print('----------------------------------------------------------------------------------')

        print()


        image.save(IMAGE_SAVE_PATH + str(file_count)+ '-' + zapis[0]+'.jpg', "JPEG")
        image.close()
        file_count = file_count + 1



    else:




        image = Image.open(IMAGE_FILE_ANKER)
        print('----------------------------Расчет для опоры '+zapis[0]+'----------------------------')
        print('Длина сегмента слева '+toFixed(zapis[2] , 2) +' м, справа '+ toFixed(zapis[3] , 2  )+' м')

        draw = ImageDraw.Draw(image)
        draw.text((1013, 177),zapis[0][6:],(0,0,0),font=font)                                                    # Печатает номер опоры
        draw.text((228, 408),toFixed(zapis[2] , 2),(0,0,0),font=font)                                            # Печатает длину пролета слева
        draw.text((592, 408),toFixed(zapis[3] , 2  ),(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((380, 525),'ДОТс-П-3кН',(0,0,0),font=font)                                                     # Печатает марку кабеля

        ## Рассчет для групп предельных состояний

        ## Идрисовская поебота 1

        G_n_L = m * g  * zapis[2] / 1000
        G_n_R = m * g  * zapis[3] / 1000

        G_nI_II_L = []
        G_nI_II_R = []
        for Y_d_i in Y_d:
            G_nI_II_L.append( toFixed((G_n_L * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
            G_nI_II_R.append( toFixed((G_n_R * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
        #print('Нормальный режим : Горизонтальная поперечная нагрузка, Н', [ 0, 0] )
        draw.text((306, 895),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 895),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1040),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1040),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1095),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1095),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1153),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1153),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1300),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1300),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1351),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1351),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        #print('Нормальный режим : Вертикальная нагрузка, Н             ', G_nI_II )


        draw.text((515, 895),G_nI_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1153),G_nI_II_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 932),G_nI_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1194),G_nI_II_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 974),G_nI_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1232),G_nI_II_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1043),G_nI_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1300),G_nI_II_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа


        draw.text((614, 895),G_nI_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1153),G_nI_II_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 932),G_nI_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1194),G_nI_II_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 974),G_nI_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1232),G_nI_II_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1043),G_nI_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1300),G_nI_II_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Нормальный режим : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
##        draw.text((513, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((925, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((513, 1050),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((925, 1050),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((513, 1145),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((925, 1145),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((513, 1220),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((925, 1220),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((513, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((925, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((319, 1312),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((319, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((701, 1312),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((701, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((608, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((1026, 950),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((608, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((1026, 1380),'0',(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('Нормальный режим : Момент силы на основание опоры, Н    ', [ 0, 0] )
        print('----------------------------------------------------------------------------------')



        fd = 0

        for yg in W_vetr_mas:
            P_vetr_ygol_L = yg * zapis[2]
            P_vetr_ygol_R = yg * zapis[3]

            P_vetr_ygol_II_L = []
            P_vetr_ygol_II_R = []
            for Y_d_i in Y_d:
                P_vetr_ygol_II_L.append( toFixed((P_vetr_ygol_L * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
                P_vetr_ygol_II_R.append( toFixed((P_vetr_ygol_R * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
            #print(toFixed(P_vetr_ygol,2))
            mom_mas = []
            if fd == 0:
                #print('Режим максимального ветра под углом 45 градусов к линии : Горизонтальная поперечная нагрузка, Н' , P_vetr_ygol_II)
                draw.text((306, 932),P_vetr_ygol_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                P_vetr_45_I_L = P_vetr_ygol_II_L[0]
                draw.text((306, 1194),P_vetr_ygol_II_L[1],(0,0,0),font=font)
                P_vetr_45_II_L = P_vetr_ygol_II_L[1]
                draw.text((411, 932),P_vetr_ygol_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                P_vetr_45_I_R = P_vetr_ygol_II_R[0]
                draw.text((411, 1194),P_vetr_ygol_II_R[1],(0,0,0),font=font)
                P_vetr_45_II_R = P_vetr_ygol_II_R[1]
##                print('Режим максимального ветра под углом 45 градусов к линии : Вертикальная нагрузка, Н             ', G_nI_II )
##                print('Режим максимального ветра под углом 45 градусов к линии : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
##                print('Режим максимального ветра под углом 45 градусов к линии : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_II) )





##                draw.text((608, 1050),moment_sili(P_vetr_ygol_II)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
##                draw.text((1026, 1050),moment_sili(P_vetr_ygol_II)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('----------------------------------------------------------------------------------')
            elif fd == 1:
##                print('Режим максимального ветра перпендикулярного линии : Горизонтальная поперечная нагрузка, Н' , P_vetr_ygol_II)

                P_vetr_90_I_L = P_vetr_ygol_II_L[0]
                P_vetr_90_II_L = P_vetr_ygol_II_L[1]
                P_vetr_90_I_R = P_vetr_ygol_II_R[0]
                P_vetr_90_II_R = P_vetr_ygol_II_R[1]



                draw.text((306, 974),P_vetr_ygol_II_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((306, 1232),P_vetr_ygol_II_L[1],(0,0,0),font=font)

                draw.text((411, 974),P_vetr_ygol_II_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
                draw.text((411, 1232),P_vetr_ygol_II_R[1],(0,0,0),font=font)
##                print('Режим максимального ветра перпендикулярного линии : Вертикальная нагрузка, Н             ', G_nI_II )
##                print('Режим максимального ветра перпендикулярного линии : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
##                print('Режим максимального ветра перпендикулярного линии : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_II) )
##                draw.text((608, 1145),moment_sili(P_vetr_ygol_II)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
##                draw.text((1026, 1145),moment_sili(P_vetr_ygol_II)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
                print('----------------------------------------------------------------------------------')
            fd = fd + 1

        ## Расчет нагрузкина опору режим гололеда с ветром
        ## вертикальная нагрузка:

        W_vertikalnaya_gololed_L =  W_gol * zapis[2]
        W_vertikalnaya_gololed_R =  W_gol * zapis[3]
        W_vertikalnaya_gololed_mas_L = []
        W_vertikalnaya_gololed_mas_R = []
        for Y_d_i in Y_d:
            W_vertikalnaya_gololed_mas_L.append( toFixed((W_vertikalnaya_gololed_L * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
            W_vertikalnaya_gololed_mas_R.append( toFixed((W_vertikalnaya_gololed_R * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1


        ## Горизонтальная нагрузка в гололед:

        P_vetr_ygol_gol_L = W_vetr_ygol * zapis[2]
        P_vetr_ygol_gol_R = W_vetr_ygol * zapis[3]
        P_vetr_ygol_gol_mas_L = []
        P_vetr_ygol_gol_mas_R = []
        for Y_d_i in Y_d:
            P_vetr_ygol_gol_mas_L.append( toFixed((P_vetr_ygol_gol_L * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
            P_vetr_ygol_gol_mas_R.append( toFixed((P_vetr_ygol_gol_R * Y_nw * Y_f_g * Y_d_i * Y_p),2) )  # Вертикальная нагрузка, Н , нормальный режим, группа 1
##        print('Режим гололеда с ветром : Горизонтальная поперечная нагрузка, Н', P_vetr_ygol_gol_mas)
        draw.text((306, 1003),P_vetr_ygol_gol_mas_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((306, 1263),P_vetr_ygol_gol_mas_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа

        draw.text((411, 1003),P_vetr_ygol_gol_mas_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((411, 1263),P_vetr_ygol_gol_mas_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
##        print('Режим гололеда с ветром : Вертикальная нагрузка, Н             ', W_vertikalnaya_gololed_mas )
        draw.text((515, 1003),W_vertikalnaya_gololed_mas_L[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((515, 1263),W_vertikalnaya_gololed_mas_L[1],(0,0,0),font=font)                                          # Печатает длину пролета справа

        draw.text((614, 1003),W_vertikalnaya_gololed_mas_R[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
        draw.text((614, 1263),W_vertikalnaya_gololed_mas_R[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
##        print('Режим гололеда с ветром : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
##        print('Режим гололеда с ветром : Момент силы на основание опоры, Н    ', moment_sili(P_vetr_ygol_gol_mas) )
##        draw.text((608, 1220),moment_sili(P_vetr_ygol_gol_mas)[0],(0,0,0),font=font)                                          # Печатает длину пролета справа
##        draw.text((1026, 1220),moment_sili(P_vetr_ygol_gol_mas)[1],(0,0,0),font=font)                                          # Печатает длину пролета справа
        print('----------------------------------------------------------------------------------')


        ## Горизонтальная продольная нагрузка при обрыве    ***** нужно вычислять исходя из максимального пролета для данного столба




        tyaga_L = []
        tyaga_R = []
        inc = 0

        for Y_d_i in Y_d:
            tyaga_L.append( toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[2] ) *  Y_d_i * Y_f_v[inc]), 2)  )
            tyaga_R.append( toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[3] ) *  Y_d_i * Y_f_v[inc]), 2)  )
            draw.text((227, 467),toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[2] )), 2),(0,0,0),font=font)
            draw.text((594, 467),toFixed((Nagruzka_ot_dlini(Provis(L_ot_T(10)), zapis[3] )), 2),(0,0,0),font=font)
            inc = inc + 1

        E1 = float(tyaga_L[0])
        F1 = float(tyaga_R[0])
        E7 = float(tyaga_L[1])
        F7 = float(tyaga_R[1])
        E4 = float(tyaga_L[0]) + 1.3 * float(W_vertikalnaya_gololed_mas_L[0])
        F4 = float(tyaga_R[0]) + 1.3 * float(W_vertikalnaya_gololed_mas_R[0])
        E10 = float(tyaga_L[1]) + 1.1 * float(W_vertikalnaya_gololed_mas_L[1])
        F10 = float(tyaga_R[1]) + 1.1 * float(W_vertikalnaya_gololed_mas_R[1])
        E2 = float(tyaga_L[0]) + 1.3 * float(P_vetr_45_I_L)
        F2 = float(tyaga_R[0]) + 1.3 * float(P_vetr_45_I_R)
        E3 = float(tyaga_L[0]) + 1.3 * float(P_vetr_90_I_L)
        F3 = float(tyaga_R[0]) + 1.3 * float(P_vetr_90_I_R)
        E8 = float(tyaga_L[1]) + 1.1 * float(P_vetr_45_II_L)
        F8 = float(tyaga_R[1]) + 1.1 * float(P_vetr_45_II_R)
        E9 = float(tyaga_L[1]) + 1.1 * float(P_vetr_90_II_L)
        F9 = float(tyaga_R[1]) + 1.1 * float(P_vetr_90_II_R)


        draw.text((711, 895),toFixed(E1),(0,0,0),font=font)
        draw.text((811, 895),toFixed(F1),(0,0,0),font=font)
        draw.text((711, 1153),toFixed(E7),(0,0,0),font=font)
        draw.text((811, 1153),toFixed(F7),(0,0,0),font=font)

        draw.text((711, 1003), toFixed(E4,2) ,(0,0,0),font=font)
        draw.text((811, 1003), toFixed(F4,2),(0,0,0),font=font)
        draw.text((711, 1263), toFixed(E10,2) ,(0,0,0),font=font)
        draw.text((811, 1263), toFixed(F10,2),(0,0,0),font=font)
        draw.text((711, 932), toFixed(E2 ,2),(0,0,0),font=font)
        draw.text((811, 932), toFixed(F2 ,2), (0,0,0),font=font)
        draw.text((711, 974), toFixed(E3 ,2),(0,0,0),font=font)
        draw.text((811, 974), toFixed(F3 ,2), (0,0,0),font=font)
        draw.text((711, 1194), toFixed(E8 ,2),(0,0,0),font=font)
        draw.text((811, 1194), toFixed(F8 ,2), (0,0,0),font=font)
        draw.text((711, 1232), toFixed(E9 ,2),(0,0,0),font=font)
        draw.text((811, 1232), toFixed(F9 ,2), (0,0,0),font=font)




##        draw.text((711, 932), str(float(tyaga_L[0]) + 1.3 * float(P_vetr_ygol_II_L[0])) ,(0,0,0),font=font)
##        draw.text((811, 932), str(float(tyaga_R[0]) + 1.3 * float(P_vetr_ygol_II_R[0])),(0,0,0),font=font)

##        draw.text((711, 932), str(float(tyaga_L[1]) + 1.1 * float(P_vetr_ygol_II_L[1])) ,(0,0,0),font=font)
##        draw.text((811, 932), str(float(tyaga_R[1]) + 1.1 * float(P_vetr_ygol_II_R[1])),(0,0,0),font=font)
##




##        print('Режим одноосного тяжения (обрыв) : Горизонтальная поперечная нагрузка, Н', [ 0, 0])
##        print('Режим одноосного тяжения (обрыв) : Вертикальная нагрузка, Н             ', G_nI_II )
##        print('Режим одноосного тяжения (обрыв) : Горизонтальная продольная нагрузка, Н', d )
##        draw.text((513, 1312),d[0],(0,0,0),font=font)
##        draw.text((925, 1312),d[1],(0,0,0),font=font)
##        print('Режим одноосного тяжения (обрыв) : Момент силы на основание опоры, Н    ', moment_sili(d) )
##        draw.text((608, 1312),moment_sili(d)[0],(0,0,0),font=font)
##        draw.text((1026, 1312),moment_sili(d)[1],(0,0,0),font=font)
        print('----------------------------------------------------------------------------------')

        d = []

        if zapis[2] >= zapis[3]:

            inc1 = 0
            for Y_d_i in Y_d:

                E5 = toFixed(float(tyaga_L[0]),2)
                F5 = 0
                E11 = toFixed(float(tyaga_L[1]),2)
                F11 = 0
                E6 = toFixed(float(tyaga_L[0]),2)
                F6 = 0
                E12 = toFixed(float(tyaga_L[1]),2)
                F12 = 0

                inc1 = inc1 + 1
        else:
            inc1 = 0
            for Y_d_i in Y_d:

                E5 = 0
                F5 = toFixed(float(tyaga_R[0]),2)
                E11 = 0
                F11 = toFixed(float(tyaga_R[1]),2)
                E6 = 0
                F6 = toFixed(float(tyaga_R[0]),2)
                E12 = 0
                F12 = toFixed(float(tyaga_R[1]),2)



                inc1 = inc1 + 1




        G1 = abs(E1-F1)
        G2 = abs(E2-F2)
        G3 = abs(E3-F3)
        G4 = abs(E4-F4)
        G5 = abs(float(E5)-float(F5))
        G6 = abs(float(E6)-float(F6))
        G7 = abs(float(E7)-float(F7))
        G8 = abs(float(E8)-float(F8))
        G9 = abs(float(E9)-float(F9))
        G10 = abs(float(E10)-float(F10))
        G11 = abs(float(E11)-float(F11))
        G12 = abs(float(E12)-float(F12))


        draw.text((711, 1043), str(E5) ,(0,0,0),font=font)
        draw.text((811, 1043), str(F5) ,(0,0,0),font=font)
        draw.text((711, 1300), str(E11) ,(0,0,0),font=font)
        draw.text((811, 1300), str(F11),(0,0,0),font=font)
        draw.text((711, 1095), str(E6) ,(0,0,0),font=font)
        draw.text((811, 1095), str(F6) ,(0,0,0),font=font)
        draw.text((711, 1351), str(E12) ,(0,0,0),font=font)
        draw.text((811, 1351), str(F12),(0,0,0),font=font)

        draw.text((919, 895), str(toFixed(G1,2)) ,(0,0,0),font=font)
        draw.text((919, 932), str(toFixed(G2,2)) ,(0,0,0),font=font)
        draw.text((919, 974), str(toFixed(G3,2)) ,(0,0,0),font=font)
        draw.text((919, 1003), str(toFixed(G4,2)),(0,0,0),font=font)
        draw.text((919, 1043), str(toFixed(G5,2)) ,(0,0,0),font=font)
        draw.text((919, 1095), str(toFixed(G6,2)) ,(0,0,0),font=font)
        draw.text((919, 1153), str(toFixed(G7,2)) ,(0,0,0),font=font)
        draw.text((919, 1194), str(toFixed(G8,2)) ,(0,0,0),font=font)
        draw.text((919, 1232), str(toFixed(G9,2)) ,(0,0,0),font=font)
        draw.text((919, 1263), str(toFixed(G10,2)),(0,0,0),font=font)
        draw.text((919, 1300), str(toFixed(G11,2)) ,(0,0,0),font=font)
        draw.text((919, 1352), str(toFixed(G12,2)) ,(0,0,0),font=font)

        H1	=	G1	*	H_z
        H2	=	G2	*	H_z
        H3	=	G3	*	H_z
        H4	=	G4	*	H_z
        H5	=	G5	*	H_z
        H6	=	G6	*	H_z
        H7	=	G7	*	H_z
        H8	=	G8	*	H_z
        H9	=	G9	*	H_z
        H10	=	G10	*	H_z
        H11	=	G11	*	H_z
        H12	=	G12	*	H_z


        draw.text((1020, 895), str(toFixed(H1,2)) ,(0,0,0),font=font)
        draw.text((1020, 932), str(toFixed(H2,2)) ,(0,0,0),font=font)
        draw.text((1020, 974), str(toFixed(H3,2)) ,(0,0,0),font=font)
        draw.text((1020, 1003), str(toFixed(H4,2)),(0,0,0),font=font)
        draw.text((1020, 1043), str(toFixed(H5,2)) ,(0,0,0),font=font)
        draw.text((1020, 1095), str(toFixed(H6,2)) ,(0,0,0),font=font)
        draw.text((1020, 1153), str(toFixed(H7,2)) ,(0,0,0),font=font)
        draw.text((1020, 1194), str(toFixed(H8,2)) ,(0,0,0),font=font)
        draw.text((1020, 1232), str(toFixed(H9,2)) ,(0,0,0),font=font)
        draw.text((1020, 1263), str(toFixed(H10,2)),(0,0,0),font=font)
        draw.text((1020, 1300), str(toFixed(H11,2)) ,(0,0,0),font=font)
        draw.text((1020, 1352), str(toFixed(H12,2)) ,(0,0,0),font=font)


        ## Вес монтажника с оснасткой


        M_montajnika = 100
        P_montajnika = M_montajnika * 9.8

        print('Монтажный режим (вес монтажника с оснасткой) : Горизонтальная поперечная нагрузка, Н', [ 0, 0])
        print('Монтажный режим (вес монтажника с оснасткой) : Вертикальная нагрузка, Н             ', [ toFixed(P_montajnika + 24.24 , 2), toFixed(P_montajnika + 12.12, 2) ] )
        draw.text((515, 1095),toFixed(P_montajnika + float(G_nI_II_L[0]) , 2),(0,0,0),font=font)
        draw.text((515, 1351),toFixed(P_montajnika + float(G_nI_II_L[1]) , 2),(0,0,0),font=font)
        draw.text((614, 1095),toFixed(P_montajnika + float(G_nI_II_R[0]) , 2),(0,0,0),font=font)
        draw.text((614, 1351),toFixed(P_montajnika + float(G_nI_II_R[1]) , 2),(0,0,0),font=font)
        print('Монтажный режим (вес монтажника с оснасткой) : Горизонтальная продольная нагрузка, Н', [ 0, 0] )
        print('Монтажный режим (вес монтажника с оснасткой) : Момент силы на основание опоры, Н    ', [ 0, 0] )
        print('----------------------------------------------------------------------------------')

        print()


        image.save(IMAGE_SAVE_PATH + str(file_count)+ '-' + zapis[0]+'.jpg', "JPEG")
        image.close()
        file_count = file_count + 1