#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <vector>

#include "vec.h"

namespace photonmap
{
    namespace utility
    {
        struct has_position_impl
        {
            template <typename R>
            static constexpr std::true_type Check(decltype(std::declval<R>().position.x, std::declval<R>().position.y,
                                                           std::declval<R>().position.z, (void)0) *);
            template <typename R>
            static constexpr std::false_type Check(...);
        };

        template <class R>
        struct has_position : decltype(has_position_impl::template Check<R>(nullptr))  // NOLINT
        {
        };

        template <class R>
        inline constexpr bool has_position_v = has_position<R>::value;  // NOLINT

        // static_assert(has_position_v<T>,
        //               "データの型は T->position->x,y,z の名前を満たし、比較可能なフィールドをもつ必要があります。");

        ///--------------------------------------------------------------///
        /// KD木の構築は下記の通り
        ///
        /// 平衡kd木
        ///     分割平面の軸はDepthによるXYZのローテーション
        ///     分割平面の点は中央値ベース
        ///--------------------------------------------------------------///
        template <typename T, std::enable_if_t<has_position_v<T>, std::nullptr_t> = nullptr>
        class KDTree
        {
        public:
            struct KDTreeNode
            {
                T *data;
                KDTreeNode *left_node;
                KDTreeNode *right_node;
                int axis;
            };

            // k-NN searchのクエリ
            struct Query
            {
                double max_distance2;   // 探索の最大半径
                size_t max_search_num;  // 最大探索点数
                Vec search_position;    // 探索中心
                Vec normal;             // 探索中心における法線
                Query(const Vec &search_position_, const Vec &normal_, const double max_distance2_,
                      const size_t max_search_num_)
                    : max_distance2(max_distance2_),
                      normal(normal_),
                      max_search_num(max_search_num_),
                      search_position(search_position_)
                {
                }
            };
            // 結果のQueueに乗せるためのデータ構造。
            struct ElementForQueue
            {
                T *point;
                double distance2;
                ElementForQueue(T *point_, const double distance2_) : point(point_), distance2(distance2_) {}
                bool operator<(const ElementForQueue &b) const { return distance2 < b.distance2; }
            };
            // KNNの結果を格納するキュー
            typedef std::priority_queue<ElementForQueue, std::vector<ElementForQueue> > ResultQueue;

            KDTree() { root_node_ = nullptr; }
            KDTree(std::vector<T> &node_values) : node_values_(node_values) {}
            ~KDTree() { DeleteNodeToLeaf(root_node_); }

            size_t Size() { return node_values_.size(); }

            void CreateKDtree() { root_node_ = BalanceKDTreeNode(node_values_.begin(), node_values_.end(), 0); }

            void AddData(const T &value) { node_values_.push_back(value); }

            auto& GetData() { return node_values_; }

            void SearchNearest(typename KDTree::ResultQueue *pqueue, typename KDTree<T>::Query &query)
            {
                SearchNodeKNN(pqueue, root_node_, query);
            }

        private:
            std::vector<T> node_values_;
            KDTreeNode *root_node_;
            void SearchNodeKNN(typename KDTree<T>::ResultQueue *pqueue, KDTreeNode *node,
                               typename KDTree<T>::Query &query)
            {
                if (node == nullptr) return;
                
                const int axis = node->axis;

                double delta;
                switch (axis)
                {
                    case 0:
                        delta = query.search_position.x - node->data->position.x;
                        break;
                    case 1:
                        delta = query.search_position.y - node->data->position.y;
                        break;
                    case 2:
                        delta = query.search_position.z - node->data->position.z;
                        break;
                }

                // 対象点<->探索中心の距離が設定半径以下　かつ　対象点<->探索中心の法線方向の距離が一定以下　という条件ならその対象点格納
                const Vec dir = node->data->position - query.search_position;
                const double distance2 = dir.length_squared();
                const double dt = dot(query.normal, dir / sqrt(distance2));
                if (distance2 < query.max_distance2 && fabs(dt) <= query.max_distance2 * 0.01)
                {
                    pqueue->push(ElementForQueue(node->data, distance2));
                    if (pqueue->size() > query.max_search_num)
                    {
                        pqueue->pop();
                        query.max_distance2 = pqueue->top().distance2;
                    }
                }
                if (delta > 0.0)
                {  // みぎ
                    SearchNodeKNN(pqueue, node->right_node, query);
                    if (delta * delta < query.max_distance2)
                    {
                        SearchNodeKNN(pqueue, node->left_node, query);
                    }
                }
                else
                {  // ひだり
                    SearchNodeKNN(pqueue, node->left_node, query);
                    if (delta * delta < query.max_distance2)
                    {
                        SearchNodeKNN(pqueue, node->right_node, query);
                    }
                }
            }

            void DeleteNodeToLeaf(KDTreeNode *node)
            {
                if (node == nullptr) return;
                DeleteNodeToLeaf(node->left_node);
                DeleteNodeToLeaf(node->right_node);
                delete node;
            }

            static auto CompareKDTreeX(const T &left, const T &right) { return left.position.x < right.position.x; }
            static auto CompareKDTreeY(const T &left, const T &right) { return left.position.y < right.position.y; }
            static auto CompareKDTreeZ(const T &left, const T &right) { return left.position.z < right.position.z; }

            KDTreeNode *BalanceKDTreeNode(typename std::vector<T>::iterator begin_it,
                                          typename std::vector<T>::iterator end_it, int depth)
            {
                if (end_it - begin_it <= 0)
                {
                    return nullptr;
                }

                const int axis = depth % 3;
                switch (axis)
                {
                    case 0:
                        std::sort(begin_it, end_it, CompareKDTreeX);
                        break;
                    case 1:
                        std::sort(begin_it, end_it, CompareKDTreeY);
                        break;
                    case 2:
                        std::sort(begin_it, end_it, CompareKDTreeZ);
                        break;
                    default:
                        break;
                }

                const int median_count = std::distance(begin_it, end_it) / 2;
                auto *node = new KDTreeNode;
                node->axis = axis;
                node->data = &(*(begin_it + median_count));

                // 再帰によるKDTree構築
                node->left_node = BalanceKDTreeNode(begin_it, begin_it + median_count, depth + 1);
                node->right_node = BalanceKDTreeNode(begin_it + median_count + 1, end_it, depth + 1);
                return node;
            }
        };
    }  // namespace utility
}  // namespace photonmap