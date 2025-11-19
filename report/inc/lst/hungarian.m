function [X_opt, f_opt] = hungarian(C, debug, maximization)

    [nRows, nCols] = size(C);
    A = C;

    % --- Приведение задачи к минимизации ---
    if maximization
        maxVal = max(C, [], "all");
        A = -(C - maxVal);
        if debug
            disp("Преобразование матрицы для задачи на максимум");
            disp(A);
        end
    end

    % --- Вычитание минимальных элементов по столбцам ---
    minCol = min(A, [], 1);
    A = A - minCol;
    if debug
        disp("Вычитание минимальных элементов по столбцам:");
        disp(A);
    end

    % --- Вычитание минимальных элементов по строкам ---
    minRow = min(A, [], 2);
    A = A - minRow;
    if debug
        disp("Вычитание минимальных элементов по строкам:");
        disp(A);
    end

    % Начальная СНН
    stars = initStars(A);
    N = size(A,1);
    quotes = zeros(N);

    if debug
        disp("Начальные 0*:");
        showDebug(A, stars, quotes, [], zeros(N));
    end

    k = nnz(stars);
    if debug
        disp("k = ");
        disp(k);
    end

    % --- Основной цикл улучшений ---
    if k < nCols
        iter = 1;
        while k < nRows
            if debug
                fprintf("Итерация %d\n", iter);
            end

            markedCols = findMarkedCols(stars);
            markedRows = zeros(size(A,2));

            while true
                [found, col, row] = findFreeZero(A, markedCols, markedRows);

                if found
                    quotes(row,col) = 1;

                    if debug
                        disp("После нахождения 0'");
                        showDebug(A, stars, quotes, markedCols, markedRows);
                    end

                    [hasStar, col2] = hasStarInRow(row, stars);

                    if hasStar
                        markedRows(row) = 1;
                        markedCols(col2) = 0;

                        if debug
                            disp("Изменили выделенные строки и столбцы");
                            showDebug(A, stars, quotes, markedCols, markedRows);
                        end
                    else
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

    % --- Построение оптимального решения ---
    X_opt = zeros(nCols);
    for i = 1:nRows
        for j = 1:nCols
            if stars(i,j)
                X_opt(i,j) = 1;
            end
        end
    end

    f_opt = sum(C(logical(stars)));
end

% ---------------------------- ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ----------------------------

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

function cols = findMarkedCols(stars)
    cols = sum(stars);
end

function [found, col, row] = findFreeZero(A, markedCols, markedRows)
    N = size(A,1);
    found = false;
    col = 0; row = 0;

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

function [hasStar, row] = hasStarInCol(col, stars)
    [hasStar, row] = hasStarInRow(col, stars');
end

function A = recalcMatrix(A, markedRows, markedCols)
    N = size(A,1);
    h = inf;

    for j = 1:N
        for i = 1:N
            if ~markedCols(j) && ~markedRows(i)
                h = min(h, A(i,j));
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
