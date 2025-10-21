% Лабораторная работа №1
% Решение задачи о назначениях венгерским методом
function mv_lab1()
    
    % Параметры для отладки и выбора типа задачи
    debug = true;       % включить вывод шагов
    maximization = true;    % задача на максимум

    % Исходная матрица стоимостей
C = [3, 5, 2, 4, 8;
    10, 10, 4, 3, 6;
    5, 6, 9, 8, 3;
    6, 2, 5, 8, 4;
    5, 4, 8, 9, 3];
    
    [nRows, nCols] = size(C);
    A = C;
    
    % Если задача на максимум: 
    % находим максимальный элемент столбца, вычитаме его
    % из каждого элемента матрицы и домножаем матрицу на −1
    if maximization == true
        maxVal = 0;
        for i = 1:nRows
            for j = 1:nCols
                if C(i,j) > maxVal
                    maxVal = C(i,j); % Находим максимальный элемент матрицы
                end
            end
        end
        for i = 1:nRows
            for j = 1:nCols
                A(i,j) = C(i,j) - maxVal; % Вычитаем максимальный элемент
                % из всех элементов матрицы
            end
        end
        A = -A; % Умножаем матрицу на -1
        if debug
            disp("Преобразование матрицы для задачи на максимум");
            disp(A);
        end
    end
    
    % Вычитание минимальных элементов
    % По столбцам
    minCol = min(A, [], 1); % Находим минимальный элемент каждого столбца
    A = A - minCol; % Вычитаем 
    if debug
        disp("Вычитание минимальных элементов по столбцам:");
        disp(A);
    end
    
    % По строкам
    minRow = min(A, [], 2); % Находим минимальный элемент каждой строки
    A = A - minRow; % Вычитаем
    if debug
        disp("Вычитание минимальных элементов по строкам:");
        disp(A);
    end
    
    % Строим начальную СНН
    stars = initStars(A);
    N = size(A, 1);
    quotes = zeros(N);
    
    if debug
        disp("Начальные 0*:");
        showDebug(A, stars, quotes, [], zeros(N));
    end
    
    % Считаем количество нулей в СНН
    k = nnz(stars); % Считаем ненулевые элеметны в м-це stars
    if debug
        disp("k = ");
        disp(k);
    end
    
    % Если k < n - нужно улучшать
    if k < nCols
        if debug, disp("k нужно улучшить."); end
        
        iter = 1;
        while k < nRows
            if debug
                fprintf("Итерация %d\n", iter);
            end
            
            markedCols = findMarkedCols(stars);
            markedRows = zeros(size(A,2));
            
            while true
                % Находим 0'
                [found, col, row] = findFreeZero(A, markedCols, markedRows);
                if found
                    % Если невыделенный нуль найден
                    quotes(row,col) = 1;
                    if debug
                        disp("После нахождения 0'");
                        showDebug(A, stars, quotes, markedCols, markedRows);
                    end
                    % Проверяем есть ли с ним в одной строке 0*
                    [hasStar, col2] = hasStarInRow(row, stars);
                    if hasStar 
                        % Если есть - снимаем выделение с этого столбца
                        % и выделяем строку с 0'
                        markedRows(row) = 1;
                        markedCols(col2) = 0;
                        if debug
                            disp("Изменили выделенные строки и столбцы");
                            showDebug(A, stars, quotes, markedCols, markedRows);
                        end
                    else
                        % Если нет - строим непродолжаемую L-цепочку
                        stars = buildLchain(stars, quotes, row, col);
                        k = k + 1;
                        markedCols = findMarkedCols(stars);
                        markedRows = zeros(size(A,2));
                        if debug
                            disp("Новые 0* после построения L-цепочки:");
                            disp(stars);
                        end
                        break;
                    end
                else
                    % Если невыделенного нуля нет
                    % Находим наим эл-т h среди невыд эл-тов
                    % Вычитаем из невыд столбцов h
                    % Добавляем к выд строкам h
                    A = recalcMatrix(A, markedRows, markedCols);
                    if debug
                        disp("Пересчитали матрицу h");
                        showDebug(A, stars, quotes, markedCols, markedRows);
                    end
                end
            end
            iter = iter + 1;
        end
    end
    
    % Оптимальное решение
    X_opt = zeros(nCols);
    for i = 1:nRows
        for j = 1:nCols
            if stars(i,j)
                X_opt(i,j) = 1;
            end
        end
    end
    disp("X_opt = ");
    disp(X_opt);
    
    f_opt = 0;
    for i = 1:nRows
        for j = 1:nCols
            if stars(i,j)
                f_opt = f_opt + C(i,j);
            end
        end
    end
    disp("f_opt = ");
    disp(f_opt);
end



% Построение начальной СНН
function stars = initStars(A)
    N = size(A,1);
    stars = zeros(N);
    for j = 1:N
        for i = 1:N
            if A(i,j) == 0 && max(stars(:,j))==0 && max(stars(i,:))==0
                stars(i,j) = 1;
                break;
            end
        end
    end
end

% Поиск отмеченных столбцов
function cols = findMarkedCols(stars)
    cols = sum(stars);
end

% Поиск неотмеченного нуля
function [found, col, row] = findFreeZero(A, markedCols, markedRows)
    N = size(A,1);
    found = false; col = 0; row = 0;
    for j = 1:N
        for i = 1:N
            if A(i,j) == 0 && ~markedCols(j) && ~markedRows(i)
                found = true;
                col = j; row = i;
                return;
            end
        end
    end
end

% Проверка наличия 0* в строке
function [hasStar, col] = hasStarInRow(row, stars)
    hasStar = false; col = 0;
    for j = 1:size(stars,2)
        if stars(row,j) == 1
            hasStar = true;
            col = j;
            return;
        end
    end
end

% Проверка наличия 0* в столбце
function [hasStar, row] = hasStarInCol(col, stars)
    [hasStar, row] = hasStarInRow(col, stars');
end

% Пересчет матрицы через минимальный элемент h
function A = recalcMatrix(A, markedRows, markedCols)
    N = size(A,1);
    h = inf;
    for j = 1:N
        for i = 1:N
            if ~markedCols(j) && ~markedRows(i)
                if A(i,j) < h
                    h = A(i,j);
                end
            end
        end
    end
    
    for j = 1:N
        for i = 1:N
            if ~markedCols(j) && ~markedRows(i)
                A(i,j) = A(i,j) - h;
            elseif markedCols(j) && markedRows(i)
                A(i,j) = A(i,j) + h;
            end
        end
    end
end

% Построение L-цепочки
function newStars = buildLchain(stars, quotes, row, col)
    newStars = stars;
    curRow = row; curCol = col;
    flag = true;
    while flag
        newStars(curRow,curCol) = 1;
        [flag, curRow] = hasStarInCol(curCol, stars);
        if flag
            newStars(curRow, curCol) = 0;
            [flag, curCol] = hasStarInRow(curRow, quotes);
        end
    end
end

% Отладочный вывод
function showDebug(A, stars, quotes, markedCols, markedRows)
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if stars(i,j)
                fprintf("%d* \t", A(i,j));
            elseif quotes(i,j)
                fprintf("%d' \t", A(i,j));
            else
                fprintf("%d \t", A(i,j));
            end
        end
        if i <= length(markedRows) && markedRows(i)
            fprintf(" +\n");
        else
            fprintf("\n");
        end
    end
    for j = 1:length(markedCols)
        if markedCols(j)
            fprintf("+ \t");
        else
            fprintf(" \t");
        end
    end
    fprintf("\n");
end
