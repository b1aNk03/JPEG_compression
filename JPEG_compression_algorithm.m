clearvars
I_in = double(imread('lena.bmp'));
% I_in = double(imread('风景.bmp'));
[m, n, dim] = size(I_in);
figure
imshow(uint8(I_in))
Get_Constant_Matrix();%生成字典和变换系数
load constant.mat
%% 初始化、色彩空间转换、4:2:2降采样
Y = zeros(m, n);
U = zeros(m, n);
V = zeros(m, n);
for i = 1 : m
    for j = 1 : n
        tmp = [I_in(i, j, 1); I_in(i, j, 2); I_in(i, j, 3)];
        Y(i, j) = T_rgb2yuv(1, :) * tmp;
        U(i, j) = T_rgb2yuv(2, :) * tmp;
        V(i, j) = T_rgb2yuv(3, :) * tmp;
    end
end
U = U + 128;
V = V + 128;
U_d = Downsampling(U);
V_d = Downsampling(V);
%% 生成Z字扫描坐标顺序
global mode_dir mapping
x = 1 : 10;
mapping = 2 .^ x - 1;
mode_dir = [ 0,  1; %→
             1, -1; %↙
             1,  0; %↓
            -1,  1];%↗
dfs_z(1, 1, 1, 1, 1); %dfs_z(x, y, dir, mode, cnt) 当前位置(x,y),下次方向dir，方向模式mode，mode = ±1, 计数器cnt，结果保存在全局变量index中
%% 编码
[Y_Q, Y_code, Y_cnt, Y_len] = Block_Encode(Y, Q_L, Dictionary_L, Dictionary_DC(1, :), 0);
[U_Q, U_code, U_cnt, U_len] = Block_Encode(U_d, Q_C, Dictionary_C, Dictionary_DC(2, :), 1);
[V_Q, V_code, V_cnt, V_len] = Block_Encode(V_d, Q_C, Dictionary_C, Dictionary_DC(2, :), 1);
CR = (m * n * 8 * 3) / (Y_len + U_len + V_len);
disp(sprintf("The Compression Rate CR = %s", num2str(CR))) %#ok<DSPS>

%% 解码
[Y_decode] = Block_Decode(Y_code, Y_cnt, Q_L, Dictionary_L, Dictionary_DC(1, :), size(Y));
[U_decode] = Block_Decode(U_code, U_cnt, Q_C, Dictionary_C, Dictionary_DC(2, :), size(U_d));
[V_decode] = Block_Decode(V_code, V_cnt, Q_C, Dictionary_C, Dictionary_DC(2, :), size(V_d));
%% 升采样，色彩空间转换
% U_reconstruct = upsampling(U_decode);
% V_reconstruct = upsampling(V_decode);
U_reconstruct = imresize(U_decode, size(U), 'bilinear');
V_reconstruct = imresize(V_decode, size(V), 'bilinear');
U_reconstruct = U_reconstruct - 128;
V_reconstruct = V_reconstruct - 128;
for i = 1 : m
    for j = 1 : n
        tmp = [Y_decode(i, j); U_reconstruct(i, j); V_reconstruct(i, j)];
        I_out(i, j, 1) = T_yuv2rgb(1, :) * tmp;
        I_out(i, j, 2) = T_yuv2rgb(2, :) * tmp;
        I_out(i, j, 3) = T_yuv2rgb(3, :) * tmp;
    end
end
figure
imshow(uint8(I_out))
MSE = sum(sum(sum((I_out - I_in) .^ 2))) / (m * n * 3);
PSNR = 10 * log10(255 ^ 2 / MSE);
disp(sprintf("The Peak Signal-to-Noise Ratio PSNR = %s", num2str(PSNR))) %#ok<DSPS>

%% 函数部分
function Img_out = Downsampling(Img_in)
    [m, n] = size(Img_in);
    Img_out = zeros(floor(m / 2), n);
    for i = 1 : floor(m / 2)
        Img_out(i, :) = Img_in(i * 2 - 1, :);%采样1 3 5 7 9 ... 行
    end
end
function dfs_z(x, y, dir, mode, cnt)
    global mode_dir index
    index(cnt, 1) = x;
    index(cnt, 2) = y;
    if x == 8 && y == 8
        return;
    end
    dir_new = dir;
    if x == 8 && y == 1
        mode = -1;
        dir = 1;
        dir_new = 4;
    end
    if judge(x, y, mode_dir(dir, 1), mode_dir(dir, 2)) %只要下次的位置在边界上就修改方向
        if mode == 1
            dir_new = mod(dir, 4) + 1;
        else
            dir_new = dir - 1;
            if dir_new == 0
                dir_new = 4;
            end
        end
    end
    dfs_z(x + mode_dir(dir, 1), y + mode_dir(dir, 2), dir_new, mode, cnt + 1);
end
function flag = judge(x, y, dx, dy)
    nx = x + dx;
    ny = y + dy;
    if nx == 1 || nx == 8 || ny == 1 || ny == 8
        flag = 1;
    else
        flag = 0;
    end
end
function [Quantization_out, code, CNT, LENGTH] = Block_Encode(Img_in, Q, Dictionary, Dictionary_DC, flag)
    %% 凑齐8*8子块，扩充区域补0
    [m, n] = size(Img_in);
    m_new = ceil(m / 8) * 8;
    n_new = ceil(n / 8) * 8;
    if flag
        Quantization_out = 128 * ones(m_new, n_new);
    else
        Quantization_out = zeros(m_new, n_new);
    end
    Quantization_out(1 : m, 1 : n) = Img_in;
    %% 8*8量化
    DC = 0;
    LENGTH = 0;
    code = {};
    CNT = 0;
    for i = 1 : 8 : m_new - 7
        for j = 1 : 8 : n_new - 7
            block = Quantization_out(i : i + 7, j : j + 7);
            block = dct2(block);
            block = floor(block ./ Q + 0.5);
            [code_ij, cnt, len] = compress(block, Dictionary, Dictionary_DC, DC);%len表示压缩后码流的长度
            for k = 1 : cnt
                code{1, CNT + k} = code_ij{1, k};
                code{2, CNT + k} = code_ij{2, k};
            end
            CNT = CNT + k;
            LENGTH = LENGTH + len;
            DC = block(1, 1);
        end
    end
end
function [code, cnt, len] = compress(block, Dictionary, Dictionary_DC, DC)
    global index mapping
    [m, ~] = size(index);
    vec = zeros(1, m);
    pre = 1;%上一个非0位置
    pos = 2;%当前位置，第一个位置是DC系数
    code = cell(2, m);
    cnt = 0;%码的长度
    len = 0;
    for i = 1 : m
        vec(i) = block(index(i, 1), index(i, 2));
    end
   %% DC系数
    DIFF = vec(1) - DC;
    SSSS = ceil(log2(abs(DIFF) + 1));
    cnt = cnt + 1;
    code{1, cnt} = Dictionary_DC(SSSS + 1);
    if DIFF <= 0
        flag = 0;
    else
        flag = 1;
    end
    if ~flag
        if DIFF ~= 0
            code{2, cnt} = dec2bin(bitxor(abs(DIFF), mapping(SSSS)), SSSS);
        else
            code{2, cnt} = dec2bin(0);
        end
    else
        code{2, cnt} = dec2bin(DIFF, SSSS);
    end
    len = len + strlength(code{1, cnt}) + length(code{2, cnt});
   %% AC系数
    while 1
        if (pos == m + 1) && (vec(pos - 1) == 0)
            cnt = cnt + 1;
            code{1, cnt} = Dictionary(17, 1);%编码为EOB
            code{2, cnt} = '';
            len = len + strlength(code{1, cnt}) + length(code{2, cnt});
            break;
        end
        if pos > m
            break;
        end
        if vec(pos) == 0
            pos = pos + 1;
            continue;
        end
        NNNN = pos - pre;%NNNN - 1对应0游程长度，NNNN是在字典中的行号
        if NNNN >= 17 %NNNN - 1 >= 16 处理中间连续0的个数大于16的情况，每16个编为ZRL，剩余少于16个的补在当前非零值前面编码
            rep = floor((NNNN - 1) / 16);
            for i = 1 : rep
                cnt = cnt + 1;
                code{1, cnt} = Dictionary(17, 2);%编码为ZRL
                code{2, cnt} = '';
                len = len + strlength(code{1, cnt}) + length(code{2, cnt});
            end
            NNNN = NNNN - 16 * rep;
        end
        pre = pos;
        if vec(pos) < 0
            SSSS = ceil(log2(-vec(pos) + 1));
            flag = 0;
        else
            SSSS = ceil(log2(vec(pos) + 1));
            flag = 1;%正数标记为1，负数标记为0，便于取反
        end
        cnt = cnt + 1;
        code{1, cnt} = Dictionary(NNNN, SSSS);
        if ~flag
            code{2, cnt} = dec2bin(bitxor(abs(vec(pos)), mapping(SSSS)), SSSS);
        else
            code{2, cnt} = dec2bin(vec(pos), SSSS);
        end
        len = len + strlength(code{1, cnt}) + length(code{2, cnt});
        pos = pos + 1;
    end
end
function Img_out = Block_Decode(code, CNT, Q, Dictionary, Dictionary_DC, sz)
    m = sz(1); m_new = ceil(m / 8) * 8;
    n = sz(2); n_new = ceil(n / 8) * 8;
    global index mapping
    Img = zeros(m_new, n_new);
    block = zeros(8, 8);
    DC = 0;
    pos = 1;%非0的位置
    pos_x = 1;
    pos_y = 1;
    for i = 1 : CNT
%         disp(sprintf("i = %s", num2str(i))) %#ok<DSPS>
%         disp(sprintf("pos_x = %s   pos_y = %s", num2str(pos_x), num2str(pos_y))) %#ok<DSPS>
       %% DC系数
        if pos == 1
            SSSS = find(Dictionary_DC == code{1, i}) - 1;
            DIFF = code{2, i};
            if judge_positive(SSSS, DIFF)%正返回1
                block(1, 1) = DC + bin2dec(DIFF);
            else
                block(1, 1) = DC - bitxor(bin2dec(DIFF), mapping(SSSS));
            end
            DC = block(1, 1);
            pos = pos + 1;
            continue;
        end
       %% AC系数
        [NNNN, SSSS] = find(Dictionary == code{1, i});
        DIFF = code{2, i};%尾码
        if NNNN == 17
            if SSSS == 1 %EOB
                pos = 1;
                Img(pos_x : pos_x + 7, pos_y : pos_y + 7) = idct2(block .* Q);
                pos_y = pos_y + 8;
                if pos_y == n_new + 1
                    pos_y = 1;
                    pos_x = pos_x + 8;
                end
                block = zeros(8, 8);
            else         %ZRL
                pos = pos + 16;
            end
            continue;
        end
        pos = pos + NNNN - 1;%有NNNN - 1个0
        if judge_positive(SSSS, DIFF)
            block(index(pos, 1), index(pos, 2)) = bin2dec(DIFF);
        else
            block(index(pos, 1), index(pos, 2)) = -bitxor(bin2dec(DIFF), mapping(SSSS));
        end
        if pos == 64
            pos = 1;
            Img(pos_x : pos_x + 7, pos_y : pos_y + 7) = idct2(block .* Q);
            pos_y = pos_y + 8;
            if pos_y == n_new + 1
                pos_y = 1;
                pos_x = pos_x + 8;
            end
            block= zeros(8, 8);
            continue;
        end
        pos = pos + 1;
    end
    Img_out = Img(1 : m, 1 : n);
    Img_out(Img_out > 255) = 255;
    Img_out(Img_out < 0) = 0;
end
function flag = judge_positive(SSSS, DIFF)
    if SSSS == 0
        flag = 1;
        return
    end
    num = bin2dec(DIFF);
    if (num < 2 ^ (SSSS - 1)) || (num > 2 ^ SSSS - 1)
        flag = 0;
    else
        flag = 1;
    end
end